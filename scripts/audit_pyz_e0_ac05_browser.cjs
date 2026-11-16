// MCL E0 / AC-05: real DOM-based browser audit. Read-only localhost GETs.
// Never starts/stops PS1/PS2, changes scientific files, or approves E0.
'use strict';
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const { performance } = require('perf_hooks');

const ROUTES = Object.freeze([
  { path: '/pyz', id: 'catalog', heading: /Молекулы PYZ/i },
  { path: '/pyz/PYZ-035', id: 'pyz-035', heading: /PYZ-035/i },
  { path: '/pyz/matrix', id: 'matrix', heading: /Матрица PYZ/i },
]);
const SIZES = Object.freeze([[1440, 900], [768, 1024], [390, 844]]);
const CHECKS = Object.freeze(['HTTP_NOT_200', 'H1_MISMATCH', 'ROOT_OVERFLOW',
  'JS_OR_HTTP_5XX_ERRORS', 'NO_KEYBOARD_FOCUS', 'SCREENSHOT_FAILED', 'VIEWPORT_MISMATCH']);

function arg(name, fallback = '') {
  const i = process.argv.indexOf(name);
  return (i >= 0 && process.argv[i + 1]) ? process.argv[i + 1] : fallback;
}
function localhostUrl(input) {
  const url = new URL(input);
  if (url.protocol !== 'http:' || !['127.0.0.1', 'localhost', '[::1]'].includes(url.hostname) ||
      (url.pathname && url.pathname !== '/') || url.search || url.hash || url.username || url.password) {
    throw new Error('Browser target must be an HTTP localhost origin without credentials/path/query');
  }
  return url.origin;
}
function reasons(item) {
  // Missing measurements are NOT_RUN, not evidence of a defective site.
  if (!item.executed) return ['NOT_RUN'];
  const r = [];
  if (item.http_status !== 200) r.push('HTTP_NOT_200');
  if (item.heading_matches !== true) r.push('H1_MISMATCH');
  if (item.layout?.horizontal_overflow === true) r.push('ROOT_OVERFLOW');
  if ((item.errors || []).length) r.push('JS_OR_HTTP_5XX_ERRORS');
  if (item.keyboard_focus !== true) r.push('NO_KEYBOARD_FOCUS');
  if (item.screenshot_saved !== true) r.push('SCREENSHOT_FAILED');
  if (item.layout?.window_inner_width !== item.width) r.push('VIEWPORT_MISMATCH');
  return r;
}
function sha(b) { return crypto.createHash('sha256').update(b).digest('hex'); }
function saveJson(file, obj) { fs.writeFileSync(file, JSON.stringify(obj, null, 2) + '\n', 'utf8'); }
function tsvValue(v) { return String(v === undefined || v === null ? '' : v).replace(/[\t\r\n]/g, ' '); }
function writeTsv(file, cols, rows) {
  fs.writeFileSync(file, cols.join('\t') + '\n' + rows.map(o => cols.map(k => tsvValue(o[k])).join('\t')).join('\n') + '\n', 'utf8');
}
function outputSummary(items) {
  const executed = items.filter(x => x.executed);
  return {
    expected_page_checks: ROUTES.length * SIZES.length,
    executed_page_checks: executed.length,
    screenshots_saved: items.filter(x => x.screenshot_saved).length,
    root_overflow_cases: executed.filter(x => x.layout?.horizontal_overflow === true).length,
    js_error_cases: executed.filter(x => (x.errors || []).length > 0).length,
    failed_technical_checks: executed.filter(x => reasons(x).length > 0).length,
    keyboard_focus_cases: executed.filter(x => x.keyboard_focus === true).length,
  };
}
function collectLayout() {
  const html = document.documentElement;
  const body = document.body;
  const viewport = window.innerWidth;
  const rootWidth = Math.max(html.scrollWidth, body?.scrollWidth || 0);
  const offenders = [];
  let inspected = 0;
  const list = document.querySelectorAll('main *, [role="main"] *');
  for (const el of list) {
    if (inspected++ > 3200 || offenders.length >= 18) break;
    if (!(el instanceof HTMLElement || el instanceof SVGElement)) continue;
    const r = el.getBoundingClientRect();
    if (!r.width || !r.height || r.right <= viewport + 3) continue;
    let node = el.parentElement;
    let scoped = false;
    while (node && node !== body) {
      const s = getComputedStyle(node);
      if (/(auto|scroll)/.test(s.overflowX) && node.scrollWidth > node.clientWidth + 2) { scoped = true; break; }
      node = node.parentElement;
    }
    if (scoped) continue;
    offenders.push({tag: el.tagName.toLowerCase(), class_name: String(el.getAttribute('class') || '').slice(0, 160),
      right: Math.round(r.right), left: Math.round(r.left), width: Math.round(r.width),
      text: (el.textContent || '').replace(/\s+/g, ' ').slice(0, 110)});
  }
  return {
    window_inner_width: window.innerWidth,
    html_client_width: html.clientWidth,
    html_scroll_width: html.scrollWidth,
    body_scroll_width: body?.scrollWidth || null,
    visual_viewport_width: window.visualViewport?.width ?? null,
    horizontal_overflow: rootWidth > html.clientWidth + 3,
    overflow_px: Math.max(0, rootWidth - html.clientWidth),
    offscreen_elements_outside_scroll_containers: offenders,
  };
}
async function launchBrowser(chromium) {
  const attempts = [];
  const candidates = [
    { channel: 'msedge' }, { channel: 'chrome' },
    ...(process.platform === 'win32' ? [
      { executablePath: path.join(process.env['PROGRAMFILES(X86)'] || 'C:\\Program Files (x86)', 'Microsoft/Edge/Application/msedge.exe') },
      { executablePath: path.join(process.env.PROGRAMFILES || 'C:\\Program Files', 'Google/Chrome/Application/chrome.exe') },
    ] : []),
  ];
  for (const choice of candidates) {
    if (choice.executablePath && !fs.existsSync(choice.executablePath)) continue;
    try {
      const browser = await chromium.launch({ ...choice, headless: true,
        args: ['--no-first-run', '--disable-extensions', '--no-default-browser-check'] });
      return { browser, selected: choice.channel || choice.executablePath };
    } catch (e) { attempts.push(`${choice.channel || choice.executablePath}: ${String(e).slice(0, 260)}`); }
  }
  throw new Error('Could not launch installed Edge/Chrome. ' + attempts.join(' | '));
}
async function execute({ origin, out, toolsRoot, expectedFile }) {
  if (fs.existsSync(out)) throw new Error(`Refusing to overwrite existing report: ${out}`);
  fs.mkdirSync(path.join(out, 'screenshots'), { recursive: true });
  const output = { contract: 'MCL-E0-AC05-1.0', timestamp_utc: new Date().toISOString(),
    origin, browser_engine: null, status: 'NOT_RUN', items: [], source_dom_comparisons: [],
    manual_visual_review_required: true, scientific_acceptance: 'NOT_ACCEPTED',
    merge_permitted: false, error: '', summary: null, evidence_files: [] };
  let browser;
  try {
    const { chromium } = require(path.join(toolsRoot, 'node_modules', 'playwright-core'));
    const launched = await launchBrowser(chromium);
    browser = launched.browser;
    output.browser_engine = launched.selected;
    for (const route of ROUTES) {
      for (const [width, height] of SIZES) {
        const item = {route: route.path, width, height, executed: false,
          http_status: null, heading: null, heading_matches: null,
          layout: null, keyboard_focus: null, errors: [], screenshot: '', screenshot_saved: false,
          screenshot_sha256: null, elapsed_domcontentloaded_ms: null, failure_reasons: []};
        const context = await browser.newContext({ viewport: {width, height}, deviceScaleFactor: 1,
          isMobile: width === 390, hasTouch: width === 390, locale: 'ru-RU', reducedMotion: 'reduce' });
        const page = await context.newPage();
        try {
          const errors = [];
          page.on('pageerror', e => errors.push('PAGEERROR: ' + String(e).slice(0, 420)));
          page.on('console', e => { if (e.type() === 'error') errors.push('CONSOLE: ' + e.text().slice(0, 420)); });
          page.on('response', res => { if (res.url().startsWith(origin) && res.status() >= 500)
            errors.push('HTTP_5XX: ' + res.status() + ' ' + res.url().slice(0, 160)); });
          const t0 = performance.now();
          const response = await page.goto(origin + route.path, {waitUntil:'domcontentloaded', timeout: 45000});
          item.executed = true;
          item.elapsed_domcontentloaded_ms = Math.round(performance.now()-t0);
          item.http_status = response ? response.status() : 0;
          try {
            await page.locator('h1').first().waitFor({ state: 'visible', timeout: 25000 });
            await page.evaluate(() => document.fonts?.ready);
            // A Next dev server may render asynchronously after DOMContentLoaded.
            await page.waitForTimeout(800);
            item.heading = (await page.locator('h1').first().textContent() || '').trim();
            item.heading_matches = route.heading.test(item.heading);
          } catch (e) { errors.push('HEADING_CHECK: ' + String(e).slice(0, 220)); item.heading_matches = false; }
          item.layout = await page.evaluate(collectLayout);
          // Three Tab presses: ensure at least one visible, keyboard-focusable control.
          let focused = false;
          for(let t=0;t<3;t++) {
            await page.keyboard.press('Tab');
            focused ||= await page.evaluate(() => {
              const el = document.activeElement;
              if (!el || el === document.body || el === document.documentElement) return false;
              const style = getComputedStyle(el);
              return style.visibility !== 'hidden' && style.display !== 'none' &&
                el.getBoundingClientRect().width > 0;
            });
          }
          item.keyboard_focus = focused;
          const screenshotName = `${route.id}-${width}x${height}.png`;
          const screenshotFile = path.join(out, 'screenshots', screenshotName);
          await page.screenshot({path:screenshotFile, fullPage:true, animations:'disabled', timeout:45000});
          item.screenshot = 'screenshots/' + screenshotName;
          item.screenshot_saved = true;
          item.screenshot_sha256 = sha(fs.readFileSync(screenshotFile));
          item.errors = errors;
        } catch(e) {
          item.errors.push('RUN: ' + String(e).slice(0, 420));
        } finally {
          item.failure_reasons = reasons(item);
          output.items.push(item);
          await context.close();
        }
      }
    }
    // Optional comparison of four/five control SMILES to previously verified E0 data.
    if (expectedFile && fs.existsSync(expectedFile)) {
      const expected = JSON.parse(fs.readFileSync(expectedFile, 'utf8'));
      for (const [id, info] of Object.entries(expected)) {
        if (!/^PYZ-\d{3}$/.test(id) || !info || typeof info.smiles !== 'string') continue;
        const record = {compound_id:id, status:'NOT_RUN', heading_matches:null, smiles_matches:null, error:''};
        const page=await browser.newPage({ viewport:{width:1440,height:900} });
        try {
          const res=await page.goto(origin + '/pyz/' + id, {waitUntil:'domcontentloaded',timeout:45000});
          await page.locator('h1').first().waitFor({state:'visible',timeout:25000});
          record.heading_matches = (await page.locator('h1').first().textContent() || '').trim() === id;
          record.smiles_matches = (await page.locator('code').first().textContent() || '').trim() === info.smiles.trim();
          record.status = res?.status()===200 && record.heading_matches && record.smiles_matches ?
            'AUTO_MATCH' : 'MISMATCH_REVIEW_REQUIRED';
        } catch(e) { record.status='CHECK_FAILED'; record.error=String(e).slice(0,350); }
        finally { await page.close(); output.source_dom_comparisons.push(record); }
      }
    }
    output.summary = outputSummary(output.items);
    const domFail = output.source_dom_comparisons.some(x=>x.status!=='AUTO_MATCH');
    output.status = output.summary.executed_page_checks === 9 &&
      output.summary.failed_technical_checks === 0 && !domFail ?
      'AUTO_PASS_PENDING_VISUAL_REVIEW' : 'AUTO_WARN_PENDING_VISUAL_REVIEW';
  } catch (e) {
    output.status = 'NOT_RUN_BROWSER_UNAVAILABLE';
    output.error = String(e).slice(0, 1200);
    output.summary = outputSummary(output.items);
  } finally {
    if (browser) await browser.close();
    const cols = ['route', 'width','height','executed','http_status','heading_matches',
      'html_scroll_width','html_client_width','overflow_px','horizontal_overflow','keyboard_focus',
      'errors','failure_reasons','screenshot','screenshot_sha256'];
    const rows = output.items.map(x=>({
      ...x, html_scroll_width:x.layout?.html_scroll_width, html_client_width:x.layout?.html_client_width,
      overflow_px:x.layout?.overflow_px, horizontal_overflow:x.layout?.horizontal_overflow,
      errors:x.errors.join(' | '), failure_reasons:x.failure_reasons.join(' | '),
    }));
    writeTsv(path.join(out,'ui_diagnostics.tsv'),cols,rows);
    writeTsv(path.join(out,'source_dom_comparisons.tsv'),
      ['compound_id','status','heading_matches','smiles_matches','error'],output.source_dom_comparisons);
    const lines = ['# MCL Э0 / AC-05 — проверка реального браузера', '',
      `- Статус: ${output.status}`,
      `- Сценариев DOM: ${output.summary.executed_page_checks}/9`,
      `- Скриншотов: ${output.summary.screenshots_saved}/9`,
      `- Случаев корневого горизонтального переполнения: ${output.summary.root_overflow_cases}`,
      `- Случаев ошибок JS/HTTP 5xx: ${output.summary.js_error_cases}`,
      `- Сценариев с предупреждениями: ${output.summary.failed_technical_checks}`,
      '', 'Ни успешная автоматическая проверка, ни скриншот не заменяют ручную оценку интерфейса.',
      'Модуль ADMET и научные файлы не изменены; Э0 остаётся НЕ ПРИНЯТЫМ.',
      '', '## Следующий шаг', '',
      'Просмотреть 9 скриншотов, проверить выделенные offscreen-элементы и подтверждённые ошибки.',
      'Исправлять CSS только при воспроизводимом root overflow или конкретном доказанном обрезании.', ''];
    fs.writeFileSync(path.join(out,'AC05_REVIEW_RU.md'),lines.join('\n'),'utf8');
    const evidence = ['ui_diagnostics.tsv','source_dom_comparisons.tsv','AC05_REVIEW_RU.md'];
    for(const x of output.items)if(x.screenshot_saved)evidence.push(x.screenshot);
    output.evidence_files = evidence.map(name=>({name, sha256:sha(fs.readFileSync(path.join(out,name)))}));
    saveJson(path.join(out,'browser_results.json'),output);
  }
  console.log(JSON.stringify({status:output.status,...output.summary,manual_visual_review_required:true},null,2));
  return output;
}
async function main(){
  const origin = localhostUrl(arg('--url','http://127.0.0.1:3100'));
  const out = path.resolve(arg('--out'));
  const toolsRoot = path.resolve(arg('--tools'));
  if (!arg('--out')||!arg('--tools')) throw Error('--out and --tools are required');
  const report = await execute({origin,out,toolsRoot,expectedFile:arg('--expected','')});
  if (report.status==='NOT_RUN_BROWSER_UNAVAILABLE') process.exitCode=2;
}
if(require.main === module) main().catch(e=>{ console.error(String(e)); process.exitCode=3; });
module.exports={ROUTES,SIZES,CHECKS,localhostUrl,reasons,outputSummary,collectLayout,execute};