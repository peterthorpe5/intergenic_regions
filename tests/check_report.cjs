/** Offline browser checks for the generated scientific dashboard. */
"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const {parseArgs} = require("node:util");
const {pathToFileURL} = require("node:url");
const {chromium} = require("playwright");
const {values} = parseArgs({options: {
  report: {type: "string"}, screenshot: {type: "string"},
  "result-json": {type: "string"},
  "browser-executable": {type: "string"},
  "chromium-module": {type: "string"},
  "with-shap-scan": {type: "boolean", default: false},
}});
assert(values.report && values.screenshot && values["result-json"],
  "Use --report, --screenshot and --result-json");

(async () => {
  const settings = {headless: true,
    executablePath: values["browser-executable"]};
  if (values["chromium-module"]) {
    const loaded = require(values["chromium-module"]);
    const runtime = loaded.default || loaded;
    settings.executablePath = await runtime.executablePath();
    settings.args = runtime.args;
  }
  const browser = await chromium.launch(settings);
  try {
    const context = await browser.newContext({offline: true,
      viewport: {width: 1440, height: 1000}});
    const page = await context.newPage();
    const errors = [];
    const network = [];
    page.on("pageerror", (error) => errors.push(error.message));
    page.on("request", (request) => {
      if (/^https?:/.test(request.url())) network.push(request.url());
    });
    await page.goto(pathToFileURL(path.resolve(values.report)).href);
    assert((await page.locator(".metric").count()) >= 6);
    assert((await page.locator("figure img").count()) >= 4);
    assert(await page.locator("figure img").evaluateAll((images) =>
      images.every((image) => image.complete && image.naturalWidth > 0)));
    if (values["with-shap-scan"]) {
      assert.equal(await page.locator("img[alt^='shap_']").count(), 6);
      assert.equal(await page.locator("img[alt^='genome_']").count(), 5);
      assert.equal(await page.locator(".figure-links a[download$='.png']").count(),
        await page.locator("figure img").count());
      assert.equal(await page.locator(".figure-links a[download$='.pdf']").count(),
        await page.locator("figure img").count());
      const pending = page.waitForEvent("download");
      await page.locator("a[download='shap_beeswarm.pdf']").click();
      const download = await pending;
      assert.equal(download.suggestedFilename(), "shap_beeswarm.pdf");
      const file = fs.readFileSync(await download.path());
      assert.equal(file.subarray(0, 4).toString(), "%PDF");
      const link = page.getByRole("link", {name: "Whole-genome scan report"});
      await link.click();
      assert.equal(await page.locator("h1").textContent(),
        "Genome-wide motif screening");
      assert.equal(await page.locator("figure img").count(), 5);
      assert(await page.locator("figure img").evaluateAll((images) =>
        images.every((image) => image.complete && image.naturalWidth > 0)));
      await page.setViewportSize({width: 390, height: 844});
      assert(await page.evaluate(() =>
        document.documentElement.scrollWidth <= window.innerWidth + 1));
      await page.setViewportSize({width: 1440, height: 1000});
      await page.goBack();
    }
    const search = page.locator("input[data-table]").first();
    const candidates = page.locator("#results-0 tbody tr");
    const initialRows = await candidates.count();
    await search.fill("positive_000");
    assert.equal(await page.locator("#results-0 tbody tr:visible").count(), 1);
    await search.fill("this_identifier_does_not_exist");
    assert.equal(await page.locator("#results-0 tbody tr:visible").count(), 0);
    await search.fill("");
    assert.equal(await page.locator("#results-0 tbody tr:visible").count(), initialRows);
    await page.locator("#results-0 button[data-sort='0']").click();
    let ranks = await candidates.evaluateAll((rows) =>
      rows.map((row) => Number(row.cells[0].textContent)));
    assert.deepEqual(ranks, [...ranks].sort((a, b) => a - b));
    await page.locator("#results-0 button[data-sort='0']").click();
    ranks = await candidates.evaluateAll((rows) =>
      rows.map((row) => Number(row.cells[0].textContent)));
    assert.deepEqual(ranks, [...ranks].sort((a, b) => b - a));
    await page.locator("#results-0 button[data-sort='0']").click();
    await page.evaluate(() => window.scrollTo(0, 0));
    await page.screenshot({path: values.screenshot});
    await page.setViewportSize({width: 390, height: 844});
    const fitsMobile = await page.evaluate(() =>
      document.documentElement.scrollWidth <= window.innerWidth + 1);
    assert(fitsMobile, "Mobile viewport must not overflow horizontally");
    assert.deepEqual(errors, []);
    assert.deepEqual(network, []);
    const result = {status: "passed", candidates: initialRows,
      embeddedPlots: await page.locator("figure img").count(),
      checks: ["offline rendering", "embedded plots", "search/reset",
        "ascending/descending numerical sorting", "mobile layout",
        "no JavaScript errors", "no external requests"]};
    if (values["with-shap-scan"]) {
      result.checks.push("six official SHAP plots", "five genome plots",
        "embedded PNG/PDF downloads", "PDF content", "linked scan report");
    }
    fs.writeFileSync(values["result-json"], JSON.stringify(result, null, 2));
    process.stdout.write(JSON.stringify(result) + "\n");
  } finally {
    await browser.close();
  }
})().catch((error) => {
  process.stderr.write(error.stack + "\n");
  process.exitCode = 1;
});
