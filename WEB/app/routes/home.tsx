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
      <div className="grid">

      </div>
      <Footer />

    </div>

  );

}
