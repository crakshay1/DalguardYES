import { createRequestHandler, RouterContextProvider } from "react-router";

const requestHandler = createRequestHandler(
  () => import("virtual:react-router/server-build"),
  import.meta.env.MODE,
);

export default {
  async fetch(request, env, ctx) {
    const loadContext = new RouterContextProvider();
    (loadContext as any).cloudflare = { env, ctx };

    return requestHandler(request, loadContext);
  },
} satisfies ExportedHandler<Env>;

