import express from "express";
import { exec } from "child_process";
import crypto from "crypto";

const app = express();
app.use(express.json());

const SECRET = process.env.WEBHOOK_SECRET || "some-secure-random-string";

app.post("/webhook", (req, res) => {
  const signature = req.headers["x-hub-signature-256"];
  if (!signature) {
    console.error("Webhook request received but missing x-hub-signature-256 header");
    return res.status(401).send("No signature");
  }

  // Calculate signature
  const hmac = crypto.createHmac("sha256", SECRET);
  // Stringify req.body or use raw body. Note: express.json() parses body, so we stringify it back.
  // To avoid formatting edge cases, stringifying req.body should be exactly matching if it was parsed from raw JSON.
  const payload = JSON.stringify(req.body);
  const digest = "sha256=" + hmac.update(payload).digest("hex");

  if (signature !== digest) {
    console.error("Signature mismatch: unauthorized webhook attempt.");
    return res.status(403).send("Invalid signature");
  }

  if (req.body.ref === "refs/heads/main") {
    console.log("Push event received for main branch. Triggering deploy...");

    // Asynchronously execute deployment to avoid webhook timeout
    exec(
      "git pull origin main && npm install && npm run build && pm2 restart beps-frontend",
      (err, stdout, stderr) => {
        if (err) {
          console.error(`Deployment failed: ${err.message}`);
          return;
        }
        console.log(`Deployment stdout:\n${stdout}`);
        if (stderr) {
          console.error(`Deployment stderr:\n${stderr}`);
        }
      }
    );

    res.status(200).send("Deployment started");
  } else {
    res.status(200).send("Not main branch. Skipping build.");
  }
});

const port = process.env.WEBHOOK_PORT || 9000;
app.listen(port, () => {
  console.log(`Git Webhook Listener running on port ${port}`);
});
