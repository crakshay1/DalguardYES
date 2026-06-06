import { useEffect, useRef, useState } from "react";
import type { Route } from "./+types/home";
import { ClientOnly } from "../components/ClientOnly";
import { ThemeToggle } from "../components/ThemeToggle";
import PillNav from "../components/PillNav";
import type { MetaFunction } from "react-router";
import { useNavigate } from "react-router";
import Footer from "../components/Footer";

export const meta: MetaFunction = () => {
  return [{ title: "DalguardYES" }, { name: "Homepage", content: "" }];
};
export function loader({ context }: Route.LoaderArgs) { }

export default function Home({ loaderData }: Route.ComponentProps) {
  const videoRef = useRef<HTMLVideoElement>(null);
  const navigate = useNavigate();
  const [sequence, setSequence] = useState("");
  const [fileName, setFileName] = useState("");

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    setFileName(file.name);
    const reader = new FileReader();
    reader.onload = (event) => {
      const text = event.target?.result as string;
      if (!text) return;

      const lines = text.split("\n");
      const sequenceLines = lines.filter(line => !line.trim().startsWith(">"));
      const parsedSequence = sequenceLines.join("").replace(/[^a-zA-Z]/g, "").toUpperCase();
      setSequence(parsedSequence);
    };
    reader.readAsText(file);
  };

  useEffect(() => {
    if (videoRef.current) {
      videoRef.current.playbackRate = 0.5;
    }
  }, []);

  const navItems = [
    { label: "Homepage", href: "/" },
    { label: "Modules", href: "/modules" },
    { label: "Browse", href: "/browse" },
    { label: "Documentation", href: "/doc" },
    { label: "Predict", href: "/predictARN" },
  ];

  return (
    <div className="bg-background text-foreground min-h-screen relative">
      <div className="relative z-30">
      </div>

      <div className="fixed top-4 right-4 flex items-center space-x-4 z-40">
        <ThemeToggle />
      </div>

      <div className="video-background">
        <img
          src="/background.png"
          className="video-background-element"
          alt="Dalguard Background"
        />


        <div className="relative z-10 flex flex-col items-center justify-center min-h-screen px-8">
          <h1 className="text-4xl font-bold text-center text-white max-w-6xl mx-auto">
            DALGUARD-YES
          </h1>
          <h1 className="text-2xl font-semibold text-center text-white mt-2">

          </h1>
          <p className="text-center text-white mt-3 text-xl max-w-2xl mx-auto">
            Submit a FASTA sequence to analyze it with our tools.
          </p>

          <div className="flex justify-center mt-10 w-full max-w-2xl">
            <div className="flex flex-col w-full border border-foreground/20 rounded-lg overflow-hidden bg-white dark:bg-zinc-900 shadow-xl">
              <textarea
                value={sequence}
                onChange={(e) => setSequence(e.target.value)}
                placeholder="Select a FASTA file or paste a FASTA sequence below."
                className="w-full h-36 px-4 py-3 bg-white dark:bg-zinc-900 text-black dark:text-white outline-none resize-none font-mono text-sm"
              />
              <div className="flex justify-between items-center px-4 py-3 bg-zinc-50 dark:bg-zinc-950 border-t border-foreground/10">
                <div className="flex items-center space-x-2">
                  <label className="cursor-pointer px-4 py-2 bg-zinc-200 hover:bg-zinc-300 dark:bg-zinc-800 dark:hover:bg-zinc-700 text-zinc-800 dark:text-zinc-100 rounded-lg text-xs font-semibold transition flex items-center space-x-1.5">
                    <span>Select a local FASTA</span>
                    <input
                      type="file"
                      accept=".fasta,.fa,.txt"
                      onChange={handleFileUpload}
                      className="hidden"
                    />
                  </label>
                  {fileName && (
                    <span className="text-xs text-foreground/75 truncate max-w-[180px] font-medium">
                      {fileName}
                    </span>
                  )}
                </div>

                <button
                  onClick={() => {
                    if (sequence.trim()) {
                      navigate(`/predictARN?sequence=${encodeURIComponent(sequence.trim().toUpperCase())}`);
                    }
                  }}
                  disabled={!sequence.trim()}
                  className="px-6 py-2 bg-primary text-primary-foreground font-semibold rounded-lg hover:opacity-90 transition disabled:opacity-50 disabled:cursor-not-allowed text-xs uppercase tracking-wider"
                >
                  Search
                </button>
              </div>
            </div>
          </div>
        </div>
      </div>



      <div className="bg-background/50 py-20 px-8 border-t border-foreground/10">
        <div className="max-w-3xl mx-auto text-center">
          <h2 className="text-3xl font-bold text-primary mb-6">
            About DalguardYES
          </h2>
          <p className="text-foreground/75 text-lg leading-relaxed">
            We do orthogonal Ribosome IG
          </p>
          <div className="flex justify-center gap-4 mt-8">
            <a
              href="/modules"
              className="px-6 py-3 bg-primary text-primary-foreground rounded-lg hover:opacity-90 transition"
            >
              DOC
            </a>
            <a
              href="/doc"
              className="px-6 py-3 border border-foreground/30 rounded-lg text-foreground hover:bg-foreground/10 transition"
            >
              RUN
            </a>
          </div>
        </div>
      </div>

      {/* Team */}
      <div className="py-16 px-8 bg-background border-t border-foreground/10">
        <h2 className="text-3xl font-bold text-center text-primary mb-6">

        </h2>

        <div className="grid grid-cols-4 md:grid-cols-8 gap-6 max-w-5xl mx-auto">

          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-red-100 text-red-700 dark:bg-red-950/40 dark:text-red-400 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              A
            </div>
            <h3 className="text-sm font-medium text-foreground">Akshay</h3>
            <p className="text-xs text-foreground/50 mt-1">PD</p>
          </div>

          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-zinc-100 text-zinc-700 dark:bg-zinc-800/40 dark:text-zinc-300 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              R
            </div>
            <h3 className="text-sm font-medium text-foreground">Rita</h3>
            <p className="text-xs text-foreground/50 mt-1"></p>
          </div>


          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-red-100 text-red-700 dark:bg-red-950/40 dark:text-red-400 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              G
            </div>
            <h3 className="text-sm font-medium text-foreground">Georgy</h3>
            <p className="text-xs text-foreground/50 mt-1">PD</p>
          </div>

          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-zinc-100 text-zinc-700 dark:bg-zinc-800/40 dark:text-zinc-300 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              S
            </div>
            <h3 className="text-sm font-medium text-foreground">Serena</h3>
            <p className="text-xs text-foreground/50 mt-1"></p>
          </div>


          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-red-100 text-red-700 dark:bg-red-950/40 dark:text-red-400 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              P
            </div>
            <h3 className="text-sm font-medium text-foreground">Paul</h3>
            <p className="text-xs text-foreground/50 mt-1">PD</p>
          </div>
          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-zinc-100 text-zinc-700 dark:bg-zinc-800/40 dark:text-zinc-300 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              T
            </div>
            <h3 className="text-sm font-medium text-foreground">Trevor</h3>
            <p className="text-xs text-foreground/50 mt-1"></p>
          </div>
          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-red-100 text-red-700 dark:bg-red-950/40 dark:text-red-400 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              T
            </div>
            <h3 className="text-sm font-medium text-foreground">Talissa</h3>
            <p className="text-xs text-foreground/50 mt-1"></p>
          </div>
          <div className="text-center">
            <div className="w-16 h-16 rounded-full bg-zinc-100 text-zinc-700 dark:bg-zinc-800/40 dark:text-zinc-300 flex items-center justify-center mx-auto mb-3 text-lg font-medium">
              E
            </div>
            <h3 className="text-sm font-medium text-foreground">Emma</h3>
            <p className="text-xs text-foreground/50 mt-1"></p>
          </div>


        </div>
      </div>
      <Footer />

    </div>

  );

}
