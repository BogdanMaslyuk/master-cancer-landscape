// MCL E0 browser audit. Uses playwright-core and an already installed Edge/Chrome.
// Never changes scientific files or terminates existing PS1/PS2 services.
const fs = require('fs');
const path = require('path');

function getArg(name, fallback) {
  const ix = process.argv.indexOf(name);
  return ix >= 0 && process.argv[ix + 1] ? process.argv[ix + 1] : fallback;
}
function saveJson(p, v) {
  fs.mkdirSync(path.dirname(p), { recursive: true });
  fs.writeFileSync(p, JSON.stringify(v, null, 2) + '\n', 'utf8');
}
function tsvValue(v) { return String(v == null ? '' : v).replace(/\t|\r|\n/g, ' '); }
async function main() {
  const url = getArg('--url', 'http://127.0.0.1:3100').replace(/\/+$/, '');
  const out = path.resolve(getArg('--out', '.'));
  const tools = path.resolve(getArg('--tools', '.'));
  const uiDir = path.join(out, 'ui');
  fs.mkdirSync(path.join(uiDir, 'screenshots'), { recursive: true });
  const result = { version: 'e0-browser-1.0', at_utc: new Date().toISOString(), url,
    status: 'NOT_RUN', browser: null, screenshots: [], control_dom_checks: [], reason: '',
    manual_review_required: true };
  const target = path.join(uiDir, 'browser_results.json');
  let browser;
  try {
    const { chromium } = require(path.join(tools, 'node_modules', 'playwright-core'));
    let launchError;
    for (const channel of ['msedge', 'chrome']) {
      try {
        browser = await chromium.launch({ channel, headless: true, args: ['--no-first-run', '--disable-extensions'] });
        result.browser = channel;
        break;
      } catch (e) { launchError = String(e); }
    }
    if (!browser) throw new Error('Neither installed Microsoft Edge nor Chrome could launch: ' + launchError);
    const expectedPath = path.join(out, 'ui_expected.json');
    const expected = fs.existsSync(expectedPath) ? JSON.parse(fs.readFileSync(expectedPath, 'utf8')) : {};
    const pages = [
      { route: '/pyz', key: 'catalog', heading: /Молекулы PYZ/i },
      { route: '/pyz/PYZ-035', key: 'pyz-035', heading: /PYZ-035/i },
      { route: '/pyz/matrix', key: 'matrix', heading: /Матрица PYZ/i },
    ];
    const sizes = [[1440, 900], [768, 1024], [390, 844]];
    for (const route of pages) {
      for (const [width, height] of sizes) {
        const context = await browser.newContext({ viewport: { width, height }, deviceScaleFactor: 1,
          locale: 'ru-RU', reducedMotion: 'reduce' });
        const page = await context.newPage();
        const errors = [];
        page.on('pageerror', e => errors.push('PAGEERROR ' + String(e)));
        page.on('console', e => { if (e.type() === 'error') errors.push('CONSOLE ' + e.text().slice(0, 500)); });
        page.on('response', resp => {
          if (resp.url().startsWith(url) && resp.status() >= 500) errors.push('HTTP ' + resp.status() + ' ' + resp.url());
        });
        const entry = { route: route.route, width, height, http_status: 0, heading: '',
          heading_matches: false, body_scroll_width: null, viewport_width: width,
          horizontal_overflow: null, keyboard_focus: false, console_errors: [],
          elapsed_domcontentloaded_ms: null, elapsed_ready_ms: null,
          screenshot: '', saved: false, status: 'NOT_RUN' };
        const start = performance.now();
        try {
          const response = await page.goto(url + route.route, { waitUntil: 'domcontentloaded', timeout: 45000 });
          entry.http_status = response ? response.status() : 0;
          entry.elapsed_domcontentloaded_ms = Math.round(performance.now() - start);
          await page.locator('h1').first().waitFor({ state: 'visible', timeout: 30000 });
          // This is a dev Next.js instance; do not block for networkidle indefinitely.
          await page.waitForTimeout(1200);
          entry.elapsed_ready_ms = Math.round(performance.now() - start);
          entry.heading = (await page.locator('h1').first().textContent() || '').trim();
          entry.heading_matches = route.heading.test(entry.heading);
          const scroll = await page.evaluate(() => ({
            scrollWidth: document.documentElement.scrollWidth,
            clientWidth: document.documentElement.clientWidth
          }));
          entry.body_scroll_width = scroll.scrollWidth;
          entry.horizontal_overflow = scroll.scrollWidth > scroll.clientWidth + 3;
          // Keyboard smoke check, not a WCAG audit.
          await page.keyboard.press('Tab');
          entry.keyboard_focus = await page.evaluate(() => {
            const el = document.activeElement;
            return !!el && el !== document.body && el !== document.documentElement;
          });
          // Restore focus to top to make screenshot consistent.
          await page.evaluate(() => document.activeElement && document.activeElement.blur());
          const name = `${route.key}-${width}x${height}.png`;
          const rel = 'screenshots/' + name;
          await page.screenshot({ path: path.join(uiDir, rel), fullPage: true, animations: 'disabled' });
          entry.screenshot = rel;
          entry.saved = true;
        } catch (e) { errors.push('RUN ' + String(e)); }
        entry.console_errors = errors;
        entry.status = (entry.saved && entry.http_status === 200 && entry.heading_matches &&
          !entry.horizontal_overflow && errors.length === 0 && entry.keyboard_focus) ? 'AUTO_PASS' : 'AUTO_WARN';
        result.screenshots.push(entry);
        await context.close();
      }
    }
    for (const [id, info] of Object.entries(expected)) {
      const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
      const record = { pyz: id, status: 'NOT_RUN', source_smiles_equal_ui: false, h1_equal_id: false, error: '' };
      try {
        const response = await page.goto(url + '/pyz/' + encodeURIComponent(id), {
          waitUntil: 'domcontentloaded', timeout: 45000 });
        await page.locator('h1').first().waitFor({ timeout: 25000 });
        record.h1_equal_id = (await page.locator('h1').first().textContent() || '').trim() === id;
        const domSmiles = (await page.locator('code').first().textContent() || '').trim();
        record.source_smiles_equal_ui = domSmiles === String(info.smiles || '').trim();
        record.status = response && response.status() === 200 && record.h1_equal_id && record.source_smiles_equal_ui
          ? 'AUTO_PASS' : 'AUTO_WARN';
      } catch (e) { record.status = 'AUTO_WARN'; record.error = String(e); }
      result.control_dom_checks.push(record);
      await page.close();
    }
    result.status = result.screenshots.length === 9 && result.screenshots.every(x => x.status === 'AUTO_PASS') &&
      result.control_dom_checks.every(x => x.status === 'AUTO_PASS') ? 'AUTO_PASS_MANUAL_REVIEW_PENDING' : 'AUTO_WARN_MANUAL_REVIEW_PENDING';
  } catch (e) {
    result.status = 'FAILED'; result.reason = String(e);
  } finally {
    if (browser) await browser.close();
    saveJson(target, result);
    const columns = ['route', 'width', 'height', 'status', 'http_status', 'heading', 'heading_matches',
      'horizontal_overflow', 'keyboard_focus', 'elapsed_domcontentloaded_ms', 'elapsed_ready_ms',
      'screenshot', 'saved', 'console_errors'];
    const lines = [columns.join('\t')];
    for (const x of result.screenshots) lines.push(columns.map(c => tsvValue(c === 'console_errors' ? x[c].join(' | ') : x[c])).join('\t'));
    fs.writeFileSync(path.join(uiDir, 'ui_checklist_automated.tsv'), lines.join('\n') + '\n', 'utf8');
  }
  console.log('E0 UI:', result.status, 'screenshots:', result.screenshots.filter(x => x.saved).length, '/ 9');
  // Do not turn a failed browser into a fake E0 acceptance.
  process.exitCode = result.status === 'FAILED' ? 2 : 0;
}
main().catch(e => { console.error(e); process.exitCode = 3; });