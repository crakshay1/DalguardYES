import { useEffect, useRef, useState } from "react";

declare global {
  interface Window {
    molstar?: any;
  }
}

interface MolstarViewerProps {
  url: string;
  format?: "pdb" | "mmcif" | "cif" | "sdf" | "gro";
  isBinary?: boolean;
  className?: string;
}

export function MolstarViewer({
  url,
  format = "pdb",
  isBinary = false,
  className = "",
}: MolstarViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [loaded, setLoaded] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    // If window.molstar is already loaded, update status immediately
    if (window.molstar) {
      setLoaded(true);
      return;
    }

    // 1. Inject Molstar CSS
    const linkId = "molstar-css";
    if (!document.getElementById(linkId)) {
      const link = document.createElement("link");
      link.id = linkId;
      link.rel = "stylesheet";
      link.href = "https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css";
      document.head.appendChild(link);
    }

    // 2. Inject Molstar JS
    const scriptId = "molstar-js";
    let script = document.getElementById(scriptId) as HTMLScriptElement | null;
    if (!script) {
      script = document.createElement("script");
      script.id = scriptId;
      script.src = "https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js";
      script.async = true;
      document.body.appendChild(script);
    }

    const handleLoad = () => {
      if (window.molstar) {
        setLoaded(true);
      } else {
        setError("Failed to initialize Molstar library from script.");
      }
    };

    const handleError = () => {
      setError("Failed to load Molstar library from CDN.");
    };

    script.addEventListener("load", handleLoad);
    script.addEventListener("error", handleError);

    return () => {
      script?.removeEventListener("load", handleLoad);
      script?.removeEventListener("error", handleError);
    };
  }, []);

  useEffect(() => {
    if (!loaded || !containerRef.current || !window.molstar) return;

    let viewerInstance: any = null;

    const initViewer = async () => {
      try {
        // Clear previous contents
        if (containerRef.current) {
          containerRef.current.innerHTML = "";
        }

        // Initialize Molstar Viewer
        viewerInstance = await window.molstar.Viewer.create(containerRef.current, {
          layoutIsExpanded: false,
          layoutShowControls: false,
          viewportShowAnimation: false,
          collapseLeftPanel: true,
        });

        // Load the structural file from the given URL
        await viewerInstance.loadStructureFromUrl(url, format, isBinary);
      } catch (err: any) {
        console.error("Molstar initialization error:", err);
        setError("Error rendering 3D structure: " + (err.message || err));
      }
    };

    initViewer();

    return () => {
      // Clean up viewer container contents
      if (containerRef.current) {
        containerRef.current.innerHTML = "";
      }
    };
  }, [loaded, url, format, isBinary]);

  if (error) {
    return (
      <div className="flex items-center justify-center bg-zinc-900 border border-red-500/30 text-red-400 text-sm p-4 rounded-lg h-[400px]">
        {error}
      </div>
    );
  }

  if (!loaded) {
    return (
      <div className="flex items-center justify-center bg-zinc-900 border border-foreground/10 text-primary/60 text-sm rounded-lg h-[400px]">
        <span className="animate-pulse">Loading structural viewer...</span>
      </div>
    );
  }

  return (
    <div
      ref={containerRef}
      className={`relative w-full h-[400px] border border-foreground/10 rounded-lg overflow-hidden bg-black ${className}`}
    />
  );
}
