import { createRequestHandler } from "@react-router/express";
import express from "express";

const app = express();

// Serve static assets from build/client
app.use(express.static("build/client"));

// Pass all other requests to React Router handler
app.all(
  "*all",
  createRequestHandler({
    build: () => import("./build/server/index.js"),
  })
);

const port = process.env.PORT || 3000;
app.listen(port, () => {
  console.log(`React Router App listening on http://localhost:${port}`);
});
