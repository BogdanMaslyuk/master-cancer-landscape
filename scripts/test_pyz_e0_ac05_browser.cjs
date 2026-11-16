'use strict';
const test=require('node:test');
const assert=require('node:assert/strict');
const { ROUTES, SIZES, localhostUrl, reasons, outputSummary } = require('./audit_pyz_e0_ac05_browser.cjs');
test('exactly nine route x viewport combinations',()=>{
 assert.deepEqual(ROUTES.map(x=>x.path),['/pyz','/pyz/PYZ-035','/pyz/matrix']);
 assert.deepEqual(SIZES, [[1440,900],[768,1024],[390,844]]);
});
test('refuse arbitrary remote domains / origins',()=>{
 assert.equal(localhostUrl('http://127.0.0.1:3100'),'http://127.0.0.1:3100');
 for(const value of ['https://127.0.0.1:3100','http://example.com','http://127.0.0.1:3100/foo','http://127.0.0.1:3100/?x=1','http://user:pass@127.0.0.1:3100'])
   assert.throws(()=>localhostUrl(value));
});
test('missing checks never become invented errors',()=>{
 const x={executed:false,width:390};
 assert.deepEqual(reasons(x),['NOT_RUN']);
});
test('real DOM overflows and JS errors are represented faithfully',()=>{
 const x={executed:true,width:390,http_status:200,heading_matches:true,
  layout:{window_inner_width:390,horizontal_overflow:true},keyboard_focus:true,
  screenshot_saved:true,errors:['CONSOLE: failure']};
 assert.deepEqual(reasons(x),['ROOT_OVERFLOW','JS_OR_HTTP_5XX_ERRORS']);
});
test('an internal scrolling table is not a root overflow',()=>{
 const x={executed:true,width:390,http_status:200,heading_matches:true,
  layout:{window_inner_width:390,horizontal_overflow:false},keyboard_focus:true,
  screenshot_saved:true,errors:[]};
 assert.deepEqual(reasons(x),[]);
});
test('keyboard and viewport mismatches are visible',()=>{
 const x={executed:true,width:390,http_status:200,heading_matches:true,
  layout:{window_inner_width:980,horizontal_overflow:false},keyboard_focus:false,
  screenshot_saved:true,errors:[]};
 assert.deepEqual(reasons(x),['NO_KEYBOARD_FOCUS','VIEWPORT_MISMATCH']);
});
test('no execution does not count as successful coverage',()=>{
 const x={executed:false,screenshot_saved:false};
 assert.deepEqual(outputSummary([x]),{
  expected_page_checks:9,executed_page_checks:0,screenshots_saved:0,
  root_overflow_cases:0,js_error_cases:0,failed_technical_checks:0,keyboard_focus_cases:0});
});