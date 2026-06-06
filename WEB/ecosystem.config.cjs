module.exports = {
  apps: [
    {
      name: "beps-frontend",
      script: "server.js",
      instances: "max",
      exec_mode: "cluster",
      env: {
        NODE_ENV: "production",
        PORT: 3000,
      },
    },
    {
      name: "git-webhook-listener",
      script: "webhook-listener.js",
      instances: 1,
      exec_mode: "fork",
      env: {
        NODE_ENV: "production",
        WEBHOOK_PORT: 9000,
        WEBHOOK_SECRET: "asr4ever", // Set this to match your Git hosting webhook secret
      },
    },
  ],
};
