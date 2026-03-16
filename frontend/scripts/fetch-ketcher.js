#!/usr/bin/env node

/**
 * Downloads Ketcher standalone from GitHub releases.
 *
 * Version is pinned in package.json under "ketcher.version".
 * Runs automatically via npm postinstall hook.
 *
 * Usage:
 *   node scripts/fetch-ketcher.js          # downloads if not present
 *   node scripts/fetch-ketcher.js --force  # re-downloads regardless
 */

const https = require("https");
const http = require("http");
const fs = require("fs");
const path = require("path");
const { execFileSync } = require("child_process");

const ROOT = path.resolve(__dirname, "..");
const STANDALONE_DIR = path.join(ROOT, "public", "standalone");
const VERSION_FILE = path.join(STANDALONE_DIR, ".ketcher-version");
const FORCE = process.argv.includes("--force");
const MAX_REDIRECTS = 5;

function getVersion() {
  const pkgPath = path.join(ROOT, "package.json");
  const pkg = JSON.parse(fs.readFileSync(pkgPath, "utf8"));
  const version = pkg.ketcher && pkg.ketcher.version;
  if (!version) {
    console.error(
      'Error: "ketcher.version" not found in package.json.\n' +
        'Add: "ketcher": { "version": "3.12.0" }',
    );
    process.exit(1);
  }
  return version;
}

function isAlreadyPresent(version) {
  if (FORCE) return false;
  if (!fs.existsSync(VERSION_FILE)) return false;
  const current = fs.readFileSync(VERSION_FILE, "utf8").trim();
  return current === version;
}

function download(url, redirectCount) {
  if (redirectCount === undefined) redirectCount = 0;
  if (redirectCount > MAX_REDIRECTS) {
    return Promise.reject(new Error("Too many redirects"));
  }

  return new Promise(function (resolve, reject) {
    var client = url.startsWith("https") ? https : http;
    client
      .get(
        url,
        { headers: { "User-Agent": "cheminformatics-microservice" } },
        function (res) {
          if (
            res.statusCode >= 300 &&
            res.statusCode < 400 &&
            res.headers.location
          ) {
            download(res.headers.location, redirectCount + 1)
              .then(resolve)
              .catch(reject);
            return;
          }
          if (res.statusCode !== 200) {
            if (res.statusCode === 404) {
              reject(
                new Error(
                  "Ketcher release not found at " +
                    url +
                    "\nCheck that the version in package.json matches a GitHub release.",
                ),
              );
            } else {
              reject(new Error("HTTP " + res.statusCode + " for " + url));
            }
            return;
          }

          var contentLength = parseInt(res.headers["content-length"], 10);
          var downloaded = 0;
          var chunks = [];
          var lastPercent = -1;

          res.on("data", function (chunk) {
            chunks.push(chunk);
            downloaded += chunk.length;
            if (contentLength) {
              var percent = Math.floor((downloaded / contentLength) * 100);
              if (percent !== lastPercent && percent % 10 === 0) {
                process.stdout.write("  " + percent + "%\r");
                lastPercent = percent;
              }
            }
          });
          res.on("end", function () {
            if (contentLength) process.stdout.write("  100%\n");
            resolve(Buffer.concat(chunks));
          });
          res.on("error", reject);
        },
      )
      .on("error", reject);
  });
}

function verify() {
  var required = [
    path.join(STANDALONE_DIR, "index.html"),
    path.join(STANDALONE_DIR, "static"),
  ];
  for (var i = 0; i < required.length; i++) {
    if (!fs.existsSync(required[i])) {
      throw new Error("Verification failed: " + required[i] + " not found");
    }
  }
}

async function main() {
  var version = getVersion();

  if (isAlreadyPresent(version)) {
    console.log(
      "Ketcher v" + version + " already present, skipping download.",
    );
    return;
  }

  var url =
    "https://github.com/epam/ketcher/releases/download/v" +
    version +
    "/ketcher-standalone-" +
    version +
    ".zip";
  console.log("Downloading Ketcher v" + version + "...");
  console.log("  URL: " + url);

  var zipBuffer = await download(url);
  console.log(
    "Downloaded " + (zipBuffer.length / 1024 / 1024).toFixed(1) + " MB",
  );

  // Clean existing standalone directory
  if (fs.existsSync(STANDALONE_DIR)) {
    fs.rmSync(STANDALONE_DIR, { recursive: true, force: true });
  }
  fs.mkdirSync(STANDALONE_DIR, { recursive: true });

  // Write zip to temp file
  var zipPath = path.join(ROOT, ".ketcher-download.zip");
  fs.writeFileSync(zipPath, zipBuffer);

  // Extract using unzip (available on macOS, Linux, CI runners)
  try {
    execFileSync("unzip", ["-qo", zipPath, "-d", STANDALONE_DIR], {
      stdio: "pipe",
    });
  } catch (err) {
    // Clean up zip on failure
    if (fs.existsSync(zipPath)) fs.unlinkSync(zipPath);
    throw new Error(
      "Extraction failed. Is 'unzip' installed?\n" +
        "  macOS: available by default\n" +
        "  Ubuntu/Debian: apt-get install unzip\n" +
        "  Alpine: apk add unzip\n" +
        "  Original error: " +
        (err.stderr ? err.stderr.toString() : err.message),
    );
  }

  // Clean up zip
  fs.unlinkSync(zipPath);

  // The zip may extract into a subdirectory — flatten if needed
  var entries = fs.readdirSync(STANDALONE_DIR);
  if (
    entries.length === 1 &&
    fs.statSync(path.join(STANDALONE_DIR, entries[0])).isDirectory()
  ) {
    var subdir = path.join(STANDALONE_DIR, entries[0]);
    var items = fs.readdirSync(subdir);
    for (var i = 0; i < items.length; i++) {
      fs.renameSync(
        path.join(subdir, items[i]),
        path.join(STANDALONE_DIR, items[i]),
      );
    }
    fs.rmdirSync(subdir);
  }

  verify();
  fs.writeFileSync(VERSION_FILE, version);
  console.log("Ketcher v" + version + " installed to public/standalone/");
}

main().catch(function (err) {
  console.error("Failed to fetch Ketcher: " + err.message);
  process.exit(1);
});
