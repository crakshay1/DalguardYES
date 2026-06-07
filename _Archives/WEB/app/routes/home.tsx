import { useEffect, useRef, useState, useMemo } from "react";
import type { Route } from "./+types/home";
import { ThemeToggle } from "../components/ThemeToggle";
import { parseFASTA, ingestJSON, ingestCSV } from "../utils/dataIngestion";
import type { MetaFunction } from "react-router";
import { useNavigate } from "react-router";
import Footer from "../components/Footer";
import {
  Trophy,
  Ruler,
  Target,
  Droplet,
  Play,
  Copy,
  Check,
  ChevronDown,
  Info,
  CheckCircle2,
  Search,
  Layers,
  Upload,
  Database,
} from "lucide-react";

export const meta: MetaFunction = () => {
  return [{ title: "RiboGuard AI - Orthogonal RBS Designer" }, { name: "Homepage", content: "Design orthogonal RBS sequences" }];
};

export function loader({ context }: Route.LoaderArgs) { }



// Initial Dashboard Datasets
const INITIAL_CANDIDATES = [
  {
    rbs: "GGAAGGA",
    spacer: "UACGUU",
    orthScore: "0.92",
    wtLeakage: "0.03",
    rbsAccess: "0.88",
    structure: "((((((.......))),). ... ((...)) AUG .....",
  },
  {
    rbs: "GGAGGAG",
    spacer: "UACUUA",
    orthScore: "0.89",
    wtLeakage: "0.04",
    rbsAccess: "0.85",
    structure: "((((((.......))))). ... ((...)) AUG .....",
  },
  {
    rbs: "AGGAGGA",
    spacer: "ACGUUA",
    orthScore: "0.86",
    wtLeakage: "0.05",
    rbsAccess: "0.82",
    structure: "(((((.........))))). ... ((...)) AUG .....",
  },
  {
    rbs: "GAGAGGA",
    spacer: "UACGCU",
    orthScore: "0.83",
    wtLeakage: "0.06",
    rbsAccess: "0.78",
    structure: "(((((.........))))). ... ((...)) AUG .....",
  },
  {
    rbs: "GGAAAGA",
    spacer: "ACGUUU",
    orthScore: "0.80",
    wtLeakage: "0.07",
    rbsAccess: "0.76",
    structure: "((((...........)))). ... ((...)) AUG .....",
  },
];

const INITIAL_SCATTER_POINTS = (() => {
  const points = [];
  const seedRandom = (s: number) => {
    const mask = 0xffffffff;
    let m_w = (123456789 + s) & mask;
    let m_z = (987654321 - s) & mask;
    return () => {
      m_z = (36969 * (m_z & 65535) + (m_z >> 16)) & mask;
      m_w = (18000 * (m_w & 65535) + (m_w >> 16)) & mask;
      let result = ((m_z << 16) + (m_w & 65535)) & mask;
      result = result / 4294967296;
      return result + 0.5;
    };
  };

  const rnd = seedRandom(42);

  for (let i = 0; i < 150; i++) {
    const logVal = -4.0 + rnd() * 4.0;
    const wtLeakage = Math.pow(10, logVal);

    const baseBinding = 0.85 - ((logVal + 4) / 4) * 0.48;
    const noise = (rnd() - 0.5) * 0.28;
    const binding = Math.max(0.05, Math.min(0.98, baseBinding + noise));
    const access = Math.max(0.02, Math.min(0.99, binding * 0.65 + rnd() * 0.45));

    points.push({
      id: i,
      wtLeakage,
      binding,
      access,
    });
  }

  for (let i = 0; i < 15; i++) {
    const logVal = -4.0 + rnd() * 1.0;
    const wtLeakage = Math.pow(10, logVal);
    const binding = 0.60 + rnd() * 0.28;
    const access = 0.75 + rnd() * 0.25;
    points.push({
      id: 1000 + i,
      wtLeakage,
      binding,
      access,
    });
  }

  return points;
})();

interface Candidate {
  rbs: string;
  spacer: string;
  orthScore: string;
  wtLeakage: string;
  rbsAccess: string;
  structure: string;
}

interface FitnessPoint {
  generation: number;
  best: number;
  avg: number;
}

interface ScatterPoint {
  id: number;
  wtLeakage: number;
  binding: number;
  access: number;
}

export default function Home({ loaderData }: Route.ComponentProps) {
  const navigate = useNavigate();

  // Input states
  const [antiSD, setAntiSD] = useState("AUCGCAAAACGAUCGUU");
  const [beforeRBS, setBeforeRBS] = useState("CUCUCUCUCUCU");
  const [afterRBS, setAfterRBS] = useState("GAGAGAGAGAGAG");
  const [wtAntiSD, setWtAntiSD] = useState("ACCUCCUUA");
  const [targetExpression, setTargetExpression] = useState("High");

  // Input file name displays
  const [antiSDFileName, setAntiSDFileName] = useState("");
  const [beforeRBSFileName, setBeforeRBSFileName] = useState("");
  const [afterRBSFileName, setAfterRBSFileName] = useState("");

  // Selection state
  const [selectedCandidateIndex, setSelectedCandidateIndex] = useState(0);

  // Copy success indicator
  const [copiedField, setCopiedField] = useState<string | null>(null);

  // Optimization run simulation
  const [isOptimizing, setIsOptimizing] = useState(false);
  const [optProgress, setOptProgress] = useState(0);

  // Stateful datasets
  const [candidates, setCandidates] = useState<Candidate[]>([]);
  const [fitnessData, setFitnessData] = useState<FitnessPoint[]>([]);
  const [scatterPoints, setScatterPoints] = useState<ScatterPoint[]>([]);

  // Parse uploaded file sequence helper
  const handleSingleFileUpload = (
    e: React.ChangeEvent<HTMLInputElement>,
    setter: (val: string) => void,
    fileNameSetter: (name: string) => void
  ) => {
    const file = e.target.files?.[0];
    if (!file) return;

    fileNameSetter(file.name);
    const reader = new FileReader();
    reader.onload = (event) => {
      const text = event.target?.result as string;
      if (!text) return;
      const parsed = parseFASTA(text, file.name);
      setter(parsed.sequence);
    };
    reader.readAsText(file);
  };

  const handleCopy = (text: string, fieldId: string) => {
    navigator.clipboard.writeText(text);
    setCopiedField(fieldId);
    setTimeout(() => setCopiedField(null), 1500);
  };

  const startOptimization = async () => {
    if (isOptimizing) return;
    setIsOptimizing(true);
    setOptProgress(10);

    try {
      setOptProgress(30);
      const response = await fetch("http://localhost:8000/api/optimize", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          antiSD: antiSD,
          beforeRBS: beforeRBS,
          afterRBS: afterRBS,
          wtAntiSD: wtAntiSD,
          targetExpression: targetExpression,
        }),
      });

      setOptProgress(70);
      if (!response.ok) {
        throw new Error(`API error: ${response.statusText}`);
      }

      const data = await response.json();

      // Update state with results from FastAPI
      if (Array.isArray(data.candidates)) {
        setCandidates(data.candidates);
      }
      if (Array.isArray(data.scatterPoints)) {
        setScatterPoints(
          data.scatterPoints.map((pt: any, idx: number) => ({
            ...pt,
            id: pt.id !== undefined ? pt.id : idx,
          }))
        );
      }
      if (Array.isArray(data.fitnessData)) {
        setFitnessData(data.fitnessData);
      }

      setSelectedCandidateIndex(0);
      setOptProgress(100);
      setTimeout(() => setIsOptimizing(false), 300);

    } catch (err) {
      console.warn("Backend server not running. Falling back to local optimization simulation.", err);

      // Fallback local simulation
      let progress = 10;
      const interval = setInterval(() => {
        progress += 15;
        if (progress >= 100) {
          clearInterval(interval);
          setOptProgress(100);
          setIsOptimizing(false);
          setSelectedCandidateIndex((prevIndex) => (prevIndex + 1) % candidates.length);
        } else {
          setOptProgress(progress);
        }
      }, 150);
    }
  };

  // Dataset Import Handler
  const handleDataImport = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      const content = event.target?.result as string;
      if (!content) return;

      if (file.name.endsWith(".json")) {
        try {
          const parsed = ingestJSON(content);

          if (parsed.inputs) {
            const parsedInputs = parsed.inputs as any;
            if (parsedInputs.antiSD) setAntiSD(parsedInputs.antiSD);
            else if (parsedInputs.orthogonalAntiSD) setAntiSD(parsedInputs.orthogonalAntiSD);

            if (parsedInputs.beforeRBS) setBeforeRBS(parsedInputs.beforeRBS);

            if (parsedInputs.afterRBS) setAfterRBS(parsedInputs.afterRBS);
            else if (parsedInputs.cdsStart) setAfterRBS(parsedInputs.cdsStart);

            if (parsedInputs.wtAntiSD) setWtAntiSD(parsedInputs.wtAntiSD);
            if (parsedInputs.targetExpression) setTargetExpression(parsedInputs.targetExpression);
          }

          if (Array.isArray(parsed.candidates)) {
            setCandidates(parsed.candidates);
            setSelectedCandidateIndex(0);
          }

          if (Array.isArray(parsed.fitnessData)) {
            setFitnessData(parsed.fitnessData);
          }

          if (Array.isArray(parsed.scatterPoints)) {
            setScatterPoints(
              parsed.scatterPoints.map((pt: any, idx: number) => ({
                ...pt,
                id: pt.id !== undefined ? pt.id : idx,
              }))
            );
          }

          alert("Dashboard data successfully loaded from JSON!");
        } catch (err) {
          alert("Error parsing JSON file. Please ensure it follows the correct schema.");
        }
      } else if (file.name.endsWith(".csv")) {
        const result = ingestCSV(content);
        if (result.type === "candidates") {
          setCandidates(result.data);
          setSelectedCandidateIndex(0);
          alert(`Loaded ${result.data.length} candidates from CSV!`);
        } else if (result.type === "scatter") {
          setScatterPoints(result.data);
          alert(`Loaded ${result.data.length} landscape data points from CSV!`);
        } else {
          alert("Could not recognize CSV structure. Please check headers.");
        }
      }
    };
    reader.readAsText(file);
  };

  const selectedCandidate = candidates[selectedCandidateIndex];

  // Auto-sync computed sequence and annotations to localStorage for Structure Inspector
  useEffect(() => {
    if (selectedCandidate) {
      const up = (selectedCandidate.five_prime_flank || beforeRBS).replace(/U/g, "T").toUpperCase();
      const rbsSeq = selectedCandidate.rbs.replace(/U/g, "T").toUpperCase();
      const spacerSeq = selectedCandidate.spacer.replace(/U/g, "T").toUpperCase();
      const down = (selectedCandidate.cds_start || afterRBS).replace(/U/g, "T").toUpperCase();

      const computedSeq = up + rbsSeq + spacerSeq + down;

      const newAnnotations = [
        { name: "Upstream (beforeRBS)", start: 0, end: up.length, direction: 1, color: "#6b7280" },
        { name: "RBS Site", start: up.length, end: up.length + rbsSeq.length, direction: 1, color: "#dc2626" },
        { name: "Spacer", start: up.length + rbsSeq.length, end: up.length + rbsSeq.length + spacerSeq.length, direction: 1, color: "#eab308" }
      ];

      if (down.length > 0) {
        newAnnotations.push({
          name: "Downstream (afterRBS)",
          start: up.length + rbsSeq.length + spacerSeq.length,
          end: computedSeq.length,
          direction: 1,
          color: "#10b981"
        });
      }

      localStorage.setItem("riboguard_computed_sequence", JSON.stringify({
        name: `RiboGuard_Candidate_${selectedCandidate.rbs}_${selectedCandidate.spacer}`,
        sequence: computedSeq,
        structure: selectedCandidate.structure,
        annotations: newAnnotations
      }));
    }
  }, [selectedCandidate, beforeRBS, afterRBS]);



  // Compute log10 bounds from actual scatter data for data-driven normalization
  const scatterBounds = useMemo(() => {
    if (scatterPoints.length === 0) return { xMin: -4, xMax: 8, yMin: -4, yMax: 8 };
    const xLogs = scatterPoints.map(p => Math.log10(Math.max(1e-12, p.wtLeakage)));
    const yLogs = scatterPoints.map(p => Math.log10(Math.max(1e-12, p.binding)));
    const xMin = Math.min(...xLogs);
    const xMax = Math.max(...xLogs);
    const yMin = Math.min(...yLogs);
    const yMax = Math.max(...yLogs);
    // Add a 5% margin
    const xPad = Math.max(0.5, (xMax - xMin) * 0.05);
    const yPad = Math.max(0.5, (yMax - yMin) * 0.05);
    return { xMin: xMin - xPad, xMax: xMax + xPad, yMin: yMin - yPad, yMax: yMax + yPad };
  }, [scatterPoints]);

  const getScatterX = (leakage: number) => {
    const logVal = Math.log10(Math.max(1e-12, leakage));
    const pct = (logVal - scatterBounds.xMin) / (scatterBounds.xMax - scatterBounds.xMin);
    return 40 + Math.max(0, Math.min(1, pct)) * 360;
  };

  const getScatterY = (binding: number) => {
    const logVal = Math.log10(Math.max(1e-12, binding));
    const pct = (logVal - scatterBounds.yMin) / (scatterBounds.yMax - scatterBounds.yMin);
    return 170 - Math.max(0, Math.min(1, pct)) * 150;
  };

  // Generate nice round log10 tick values within a range
  const getLogTicks = (minLog: number, maxLog: number, maxTicks: number = 5) => {
    const ticks: number[] = [];
    const step = Math.max(1, Math.round((maxLog - minLog) / maxTicks));
    const start = Math.ceil(minLog);
    for (let e = start; e <= Math.floor(maxLog); e += step) {
      ticks.push(e);
    }
    return ticks;
  };

  const getScatterColor = (access: number) => {
    const hue = 270 - access * 210;
    return `hsl(${hue}, 80%, 50%)`;
  };

  // Interactive scatter point tooltip
  const [hoveredScatterPoint, setHoveredScatterPoint] = useState<{
    id: number;
    wtLeakage: number;
    binding: number;
    access: number;
    x: number;
    y: number;
  } | null>(null);

  const getFitnessX = (gen: number) => 40 + (gen / 24) * 440;
  const getFitnessY = (fit: number) => 170 - fit * 140;

  const bestD = useMemo(() => {
    if (fitnessData.length === 0) return "";
    const len = fitnessData.length;
    const getX = (gen: number) => 40 + (gen / (len - 1)) * 440;
    return "M " + fitnessData.map((d, idx) => `${getX(idx)},${getFitnessY(d.best)}`).join(" L ");
  }, [fitnessData]);

  const avgD = useMemo(() => {
    if (fitnessData.length === 0) return "";
    const len = fitnessData.length;
    const getX = (gen: number) => 40 + (gen / (len - 1)) * 440;
    return "M " + fitnessData.map((d, idx) => `${getX(idx)},${getFitnessY(d.avg)}`).join(" L ");
  }, [fitnessData]);

  // GA Fitness Interactive Tooltip State
  const [hoveredFitnessData, setHoveredFitnessData] = useState<{
    generation: number;
    best: number;
    avg: number;
    x: number;
    y: number;
  } | null>(null);

  const handleFitnessMouseMove = (e: React.MouseEvent<SVGSVGElement, MouseEvent>) => {
    if (fitnessData.length === 0) return;
    const rect = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const plotWidth = 440;
    const startX = 40;
    const pct = Math.max(0, Math.min(1, (x - startX) / plotWidth));
    const idx = Math.round(pct * (fitnessData.length - 1));

    const dataPoint = fitnessData[idx];
    if (dataPoint) {
      const getX = (gen: number) => 40 + (gen / (fitnessData.length - 1)) * 440;
      setHoveredFitnessData({
        ...dataPoint,
        x: getX(idx),
        y: getFitnessY(dataPoint.best),
      });
    }
  };

  const handleFitnessMouseLeave = () => {
    setHoveredFitnessData(null);
  };

  return (
    <div className="bg-background text-foreground min-h-screen flex font-sans overflow-hidden">

      {/* Sidebar - Inherits from global theme */}
      <aside className="w-64 bg-sidebar text-sidebar-foreground flex flex-col justify-between shrink-0 border-r border-sidebar-border select-none h-screen">
        <div>
          {/* Logo Header */}
          <div className="p-6 flex items-center space-x-3 border-b border-sidebar-border/60">
            <div className="bg-sidebar-primary p-1.5  text-sidebar-primary-foreground">
              <svg className="w-6 h-6" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5">
                <path d="M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10z" />
                <path d="M8 11h8M8 15h8M8 7h8" strokeLinecap="round" />
              </svg>
            </div>
            <span className="font-extrabold text-sidebar-foreground text-lg tracking-wide">RiboGuard AI</span>
          </div>

          {/* Navigation Links */}
          <nav className="p-4 space-y-1">
            {[
              { name: "Design", icon: Layers, path: "/", active: true },
              { name: "Structure Inspector", icon: Search, path: "/inspector" },
            ].map((item) => (
              <button
                key={item.name}
                onClick={() => !item.disabled && item.path && navigate(item.path)}
                className={`w-full flex items-center space-x-3 px-4 py-3  text-sm font-semibold transition-all duration-150 ${item.active
                  ? "bg-sidebar-primary text-sidebar-primary-foreground "
                  : "text-sidebar-foreground/70 hover:bg-sidebar-accent hover:text-sidebar-accent-foreground"
                  }`}
              >
                <item.icon className="w-4.5 h-4.5 shrink-0" />
                <span>{item.name}</span>
              </button>
            ))}
          </nav>
        </div>

      </aside>

      {/* Main Container */}
      <div className="flex-1 flex flex-col h-screen overflow-hidden">

        {/* Top Header */}
        <header className="px-8 py-5 border-b border-border flex items-center justify-between bg-card text-card-foreground backdrop-blur-md shrink-0">
          <div>
            <h1 className="text-2xl font-bold tracking-tight leading-tight">
              Orthogonal RBS Designer
            </h1>
            <p className="text-xs text-muted-foreground mt-1 font-medium">
              Design RBS + spacer sequences for a custom anti-Shine-Dalgarno.
            </p>
          </div>

          <div className="flex items-center space-x-4">
            {/* Theme Toggle */}



            {/* Import Data Trigger */}
            <label className="flex items-center space-x-1.5 text-xs text-foreground/80 hover:text-foreground font-bold cursor-pointer transition-colors bg-card hover:bg-muted/30 px-3.5 py-2  border border-border ">
              <Database className="w-4 h-4 text-primary" />
              <span>Import Dataset</span>
              <input
                type="file"
                accept=".json,.csv"
                className="hidden"
                onChange={handleDataImport}
              />
            </label>

            {/* Run Button */}
            <button
              onClick={startOptimization}
              disabled={isOptimizing}
              className={`flex items-center space-x-2 px-5 py-2.5  text-primary-foreground font-bold text-sm shadow-md transition-all duration-200 ${isOptimizing
                ? "bg-primary/60 cursor-not-allowed"
                : "bg-primary hover:bg-primary/90 active:scale-[0.98] shadow-primary/10"
                }`}
            >
              {isOptimizing ? (
                <>
                  <div className="w-4 h-4 border-2 border-primary-foreground border-t-transparent rounded-full animate-spin" />
                  <span>Running ({optProgress}%)</span>
                </>
              ) : (
                <>
                  <Play className="w-4 h-4 fill-current" />
                  <span>Run Optimization</span>
                </>
              )}
            </button>
          </div>
        </header>

        {/* Dashboard Canvas Area */}
        <div className="flex-1 overflow-y-auto p-8 space-y-6 bg-muted/20">
          <div className="grid grid-cols-6 gap-6 max-w-[1400px] mx-auto">

            {/* CARD A: Designer Inputs */}
            <div className="col-span-3 bg-card text-card-foreground p-6 border border-border flex flex-col justify-between space-y-4">
              <div className="border-b border-border pb-2">
                <h2 className="font-bold text-card-foreground/80 text-xs tracking-wider uppercase">
                  Designer Inputs
                </h2>
              </div>

              <div className="space-y-4">
                {/* Input 1: antiSD */}
                <div className="space-y-1.5">
                  <div className="flex justify-between items-center">
                    <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                      antiSD Sequence (--sequence)
                    </label>
                    <label className="flex items-center space-x-1.5 text-xs text-primary hover:text-primary/80 font-bold cursor-pointer transition-colors bg-primary/5 px-2.5 py-0.5 border border-primary/10">
                      <Upload className="w-3.5 h-3.5" />
                      <span>{antiSDFileName ? antiSDFileName : "Upload file"}</span>
                      <input
                        type="file"
                        accept=".fasta,.txt,.seq,.fa"
                        className="hidden"
                        onChange={(e) => handleSingleFileUpload(e, setAntiSD, setAntiSDFileName)}
                      />
                    </label>
                  </div>
                  <div className="relative flex items-center">
                    <input
                      type="text"
                      value={antiSD}
                      onChange={(e) => setAntiSD(e.target.value.toUpperCase())}
                      className="w-full bg-muted/40 border border-border px-3 py-2 text-sm font-mono focus:outline-none focus:ring-1 focus:ring-primary pr-10 uppercase text-foreground"
                      placeholder="Paste antiSD sequence or local file path"
                    />
                    <button
                      onClick={() => handleCopy(antiSD, "antiSD")}
                      className="absolute right-2 text-muted-foreground hover:text-foreground p-1"
                    >
                      {copiedField === "antiSD" ? (
                        <Check className="w-4 h-4 text-primary animate-in fade-in zoom-in-50 duration-200" />
                      ) : (
                        <Copy className="w-4 h-4" />
                      )}
                    </button>
                  </div>
                </div>

                {/* Input 2: beforeRBS */}
                <div className="space-y-1.5">
                  <div className="flex justify-between items-center">
                    <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                      beforeRBS Sequence (--mrna5)
                    </label>
                    <label className="flex items-center space-x-1.5 text-xs text-primary hover:text-primary/80 font-bold cursor-pointer transition-colors bg-primary/5 px-2.5 py-0.5 border border-primary/10">
                      <Upload className="w-3.5 h-3.5" />
                      <span>{beforeRBSFileName ? beforeRBSFileName : "Upload file"}</span>
                      <input
                        type="file"
                        accept=".fasta,.txt,.seq,.fa"
                        className="hidden"
                        onChange={(e) => handleSingleFileUpload(e, setBeforeRBS, setBeforeRBSFileName)}
                      />
                    </label>
                  </div>
                  <div className="relative flex items-center">
                    <input
                      type="text"
                      value={beforeRBS}
                      onChange={(e) => setBeforeRBS(e.target.value.toUpperCase())}
                      className="w-full bg-muted/40 border border-border px-3 py-2 text-sm font-mono focus:outline-none focus:ring-1 focus:ring-primary pr-10 uppercase text-foreground"
                      placeholder="Paste beforeRBS sequence or local file path"
                    />
                    <button
                      onClick={() => handleCopy(beforeRBS, "beforeRBS")}
                      className="absolute right-2 text-muted-foreground hover:text-foreground p-1"
                    >
                      {copiedField === "beforeRBS" ? (
                        <Check className="w-4 h-4 text-primary animate-in fade-in zoom-in-50 duration-200" />
                      ) : (
                        <Copy className="w-4 h-4" />
                      )}
                    </button>
                  </div>
                </div>

                {/* Input 3: afterRBS */}
                <div className="space-y-1.5">
                  <div className="flex justify-between items-center">
                    <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                      afterRBS Sequence (--mrna3)
                    </label>
                    <label className="flex items-center space-x-1.5 text-xs text-primary hover:text-primary/80 font-bold cursor-pointer transition-colors bg-primary/5 px-2.5 py-0.5 border border-primary/10">
                      <Upload className="w-3.5 h-3.5" />
                      <span>{afterRBSFileName ? afterRBSFileName : "Upload file"}</span>
                      <input
                        type="file"
                        accept=".fasta,.txt,.seq,.fa"
                        className="hidden"
                        onChange={(e) => handleSingleFileUpload(e, setAfterRBS, setAfterRBSFileName)}
                      />
                    </label>
                  </div>
                  <div className="relative flex items-center">
                    <input
                      type="text"
                      value={afterRBS}
                      onChange={(e) => setAfterRBS(e.target.value.toUpperCase())}
                      className="w-full bg-muted/40 border border-border px-3 py-2 text-sm font-mono focus:outline-none focus:ring-1 focus:ring-primary pr-10 uppercase text-foreground"
                      placeholder="Paste afterRBS sequence or local file path"
                    />
                    <button
                      onClick={() => handleCopy(afterRBS, "afterRBS")}
                      className="absolute right-2 text-muted-foreground hover:text-foreground p-1"
                    >
                      {copiedField === "afterRBS" ? (
                        <Check className="w-4 h-4 text-primary animate-in fade-in zoom-in-50 duration-200" />
                      ) : (
                        <Copy className="w-4 h-4" />
                      )}
                    </button>
                  </div>
                </div>

                {/* Additional parameters: WT anti-SD and Target Expression */}
                <div className="grid grid-cols-2 gap-4 pt-2">
                  <div className="space-y-1">
                    <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                      WT anti-SD
                    </label>
                    <div className="relative flex items-center">
                      <input
                        type="text"
                        value={wtAntiSD}
                        onChange={(e) => setWtAntiSD(e.target.value.toUpperCase())}
                        className="w-full bg-muted/40 border border-border px-3 py-2 text-sm font-mono pr-10 focus:outline-none focus:ring-1 focus:ring-primary uppercase text-foreground"
                      />
                      <button
                        onClick={() => handleCopy(wtAntiSD, "wt")}
                        className="absolute right-2 text-muted-foreground hover:text-foreground p-1"
                      >
                        {copiedField === "wt" ? (
                          <Check className="w-4 h-4 text-primary animate-in fade-in zoom-in-50 duration-200" />
                        ) : (
                          <Copy className="w-4 h-4" />
                        )}
                      </button>
                    </div>
                  </div>

                  <div className="space-y-1">
                    <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                      Target expression
                    </label>
                    <div className="relative">
                      <select
                        value={targetExpression}
                        onChange={(e) => setTargetExpression(e.target.value)}
                        className="w-full bg-muted/40 border border-border px-3 py-2 text-sm appearance-none focus:outline-none focus:ring-1 focus:ring-primary cursor-pointer text-foreground"
                      >
                        <option value="High">High</option>
                        <option value="Medium">Medium</option>
                        <option value="Low">Low</option>
                      </select>
                      <ChevronDown className="w-4 h-4 text-muted-foreground absolute right-3 top-3.5 pointer-events-none" />
                    </div>
                  </div>
                </div>
              </div>
            </div>

            {/* CARD B: KPI Metrics */}
            <div className="col-span-3 grid grid-cols-4 gap-4">

              {/* Best RBS */}
              <div className="bg-card text-card-foreground p-4  border border-border  flex flex-col justify-between items-center text-center">
                <div className="flex items-center space-x-1.5 text-muted-foreground">
                  <Trophy className="w-4.5 h-4.5 text-primary" />
                  <span className="text-xs font-bold uppercase tracking-wider">Best RBS</span>
                </div>
                <div className="text-primary font-extrabold text-lg my-1.5 select-all font-mono">
                  {selectedCandidate ? selectedCandidate.rbs : "N/A"}
                </div>
                <div className="text-[10px] text-muted-foreground font-semibold">
                  (Spacer {selectedCandidate ? selectedCandidate.spacer.length : 0} nt)
                </div>
              </div>

              {/* Best Spacer */}
              <div className="bg-card text-card-foreground p-4  border border-border  flex flex-col justify-between items-center text-center">
                <div className="flex items-center space-x-1.5 text-muted-foreground">
                  <Ruler className="w-4.5 h-4.5 text-primary/80" />
                  <span className="text-xs font-bold uppercase tracking-wider">Best Spacer</span>
                </div>
                <div className="text-foreground font-extrabold text-lg my-1.5">
                  {selectedCandidate ? selectedCandidate.spacer.length : 0} nt
                </div>
                <div className="text-[10px] text-muted-foreground font-semibold">
                  (Range: 3-12 nt)
                </div>
              </div>

              {/* Orth Score */}
              <div className="bg-card text-card-foreground p-4  border border-border  flex flex-col justify-between items-center text-center">
                <div className="flex items-center space-x-1.5 text-muted-foreground">
                  <Target className="w-4.5 h-4.5 text-primary" />
                  <span className="text-xs font-bold uppercase tracking-wider">Orth Score</span>
                </div>
                <div className="text-primary font-extrabold text-lg my-1.5">
                  {selectedCandidate ? selectedCandidate.orthScore : "N/A"}
                </div>
                <div className="text-[10px] text-muted-foreground font-semibold">
                  (Higher is better)
                </div>
              </div>

              {/* Native Leak */}
              <div className="bg-card text-card-foreground p-4  border border-border  flex flex-col justify-between items-center text-center">
                <div className="flex items-center space-x-1.5 text-muted-foreground">
                  <Droplet className="w-4.5 h-4.5 text-primary/80" />
                  <span className="text-xs font-bold uppercase tracking-wider">Native Leak</span>
                </div>
                <div className="text-foreground/80 font-extrabold text-lg my-1.5">
                  {selectedCandidate ? selectedCandidate.wtLeakage : "N/A"}
                </div>
                <div className="text-[10px] text-muted-foreground font-semibold">
                  (Lower is better)
                </div>
              </div>
            </div>

            {/* CARD C: GA Fitness Evolution */}
            <div className="col-span-3 bg-card text-card-foreground p-6  border border-border  relative group">
              <div className="flex items-center space-x-1.5 mb-2">
                <h3 className="font-bold text-card-foreground text-xs tracking-wider uppercase">
                  GA Fitness Evolution
                </h3>
                <div className="relative cursor-pointer group/info">
                  <Info className="w-4 h-4 text-muted-foreground hover:text-foreground" />
                  <div className="absolute left-1/2 -translate-x-1/2 bottom-full mb-2 w-48 bg-neutral-900 text-white text-[10px] p-2.5  opacity-0 pointer-events-none group-hover/info:opacity-100 transition-opacity z-50 shadow-lg leading-relaxed">
                    Tracks the maximum (Best) and mean (Average) fitness scores of RBS designs over the generations of the Genetic Algorithm.
                  </div>
                </div>
              </div>

              {/* Chart Legend */}
              <div className="flex items-center space-x-4 text-xs font-bold mb-4">
                <div className="flex items-center space-x-1.5">
                  <span className="w-3.5 h-0.5 bg-primary inline-block" />
                  <span className="text-muted-foreground">Best Fitness</span>
                </div>
                <div className="flex items-center space-x-1.5">
                  <span className="w-3.5 h-0.5 border-t border-dashed border-muted-foreground inline-block" />
                  <span className="text-muted-foreground">Average Fitness</span>
                </div>
              </div>

              {/* SVG Plot */}
              <div className="relative">
                <svg
                  viewBox="0 0 500 200"
                  className="w-full overflow-visible"
                  onMouseMove={handleFitnessMouseMove}
                  onMouseLeave={handleFitnessMouseLeave}
                >
                  {/* Grid Lines */}
                  {[0, 0.2, 0.4, 0.6, 0.8, 1.0].map((val) => (
                    <line
                      key={val}
                      x1="40"
                      y1={getFitnessY(val)}
                      x2="480"
                      y2={getFitnessY(val)}
                      className="stroke-border"
                      strokeWidth="1"
                    />
                  ))}
                  {fitnessData.length > 0 && [0, 0.25, 0.5, 0.75, 1.0].map((pct) => {
                    const idx = Math.round(pct * (fitnessData.length - 1));
                    const len = fitnessData.length;
                    const x = 40 + (idx / (len - 1)) * 440;
                    return (
                      <line
                        key={pct}
                        x1={x}
                        y1="20"
                        x2={x}
                        y2="170"
                        className="stroke-border"
                        strokeWidth="1"
                      />
                    );
                  })}

                  {/* Axes */}
                  <line x1="40" y1="170" x2="480" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />
                  <line x1="40" y1="20" x2="40" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />

                  {/* Axis Labels */}
                  {fitnessData.length > 0 && [0, 0.25, 0.5, 0.75, 1.0].map((pct) => {
                    const idx = Math.round(pct * (fitnessData.length - 1));
                    const gen = fitnessData[idx]?.generation ?? idx;
                    const len = fitnessData.length;
                    const x = 40 + (idx / (len - 1)) * 440;
                    return (
                      <text
                        key={pct}
                        x={x}
                        y="185"
                        textAnchor="middle"
                        className="text-[9px] fill-muted-foreground font-bold"
                      >
                        {gen}
                      </text>
                    );
                  })}
                  <text x="260" y="198" textAnchor="middle" className="text-[10px] fill-muted-foreground font-bold">
                    Generation
                  </text>

                  {[0, 0.2, 0.4, 0.6, 0.8, 1.0].map((val) => (
                    <text
                      key={val}
                      x="30"
                      y={getFitnessY(val) + 3}
                      textAnchor="end"
                      className="text-[9px] fill-muted-foreground font-bold"
                    >
                      {val.toFixed(1)}
                    </text>
                  ))}
                  <text
                    x="12"
                    y="95"
                    textAnchor="middle"
                    transform="rotate(-90, 12, 95)"
                    className="text-[10px] fill-muted-foreground font-bold"
                  >
                    Fitness
                  </text>

                  {/* Line Paths */}
                  {fitnessData.length > 0 && (
                    <>
                      <path
                        d={avgD}
                        fill="none"
                        className="stroke-muted-foreground transition-all duration-300"
                        strokeWidth="1.5"
                        strokeDasharray="4,3"
                      />
                      <path
                        d={bestD}
                        fill="none"
                        className="stroke-primary transition-all duration-300"
                        strokeWidth="2"
                      />
                    </>
                  )}

                  {/* Hover Elements */}
                  {hoveredFitnessData && (
                    <>
                      <line
                        x1={hoveredFitnessData.x}
                        y1="20"
                        x2={hoveredFitnessData.x}
                        y2="170"
                        className="stroke-muted-foreground"
                        strokeWidth="1"
                        strokeDasharray="2,2"
                      />
                      <circle
                        cx={hoveredFitnessData.x}
                        cy={getFitnessY(hoveredFitnessData.best)}
                        r="4.5"
                        className="fill-primary stroke-card"
                        strokeWidth="1.5"
                      />
                      <circle
                        cx={hoveredFitnessData.x}
                        cy={getFitnessY(hoveredFitnessData.avg)}
                        r="4.5"
                        className="fill-muted-foreground stroke-card"
                        strokeWidth="1.5"
                      />
                    </>
                  )}
                </svg>

                {/* Fitness Tooltip overlay */}
                {hoveredFitnessData && (
                  <div
                    className="absolute bg-neutral-900/95 text-white text-[10px] p-2.5  shadow-lg pointer-events-none z-30 flex flex-col space-y-1 border border-neutral-800"
                    style={{
                      left: `${(hoveredFitnessData.x / 500) * 100}%`,
                      top: "24px",
                      transform: hoveredFitnessData.x > 250 ? "translateX(-110%)" : "translateX(10%)",
                    }}
                  >
                    <div className="font-bold text-neutral-300 border-b border-neutral-800 pb-1 mb-1">
                      Gen {hoveredFitnessData.generation}
                    </div>
                    <div className="flex justify-between space-x-4">
                      <span className="text-neutral-400 font-semibold">Best:</span>
                      <span className="font-bold text-red-400">{hoveredFitnessData.best}</span>
                    </div>
                    <div className="flex justify-between space-x-4">
                      <span className="text-neutral-400 font-semibold">Average:</span>
                      <span className="font-bold text-neutral-400">{hoveredFitnessData.avg}</span>
                    </div>
                  </div>
                )}
              </div>
            </div>

            {/* CARD D: Orthogonality Landscape */}
            <div className="col-span-3 bg-card text-card-foreground p-6  border border-border  relative group">
              <div className="flex items-center space-x-1.5 mb-4">
                <h3 className="font-bold text-card-foreground text-xs tracking-wider uppercase">
                  Orthogonality Landscape
                </h3>
                <div className="relative cursor-pointer group/info">
                  <Info className="w-4 h-4 text-muted-foreground hover:text-foreground" />
                  <div className="absolute left-1/2 -translate-x-1/2 bottom-full mb-2 w-48 bg-neutral-900 text-white text-[10px] p-2.5  opacity-0 pointer-events-none group-hover/info:opacity-100 transition-opacity z-50 shadow-lg leading-relaxed">
                    Scatter plot representing candidate designs. Y-axis is orthogonal binding strength. X-axis is wild-type leakage (log scale). Hover over points to view metrics.
                  </div>
                </div>
              </div>

              <div className="relative">
                <svg viewBox="0 0 500 200" className="w-full overflow-visible">
                  <defs>
                    <linearGradient id="access-grad" x1="0%" y1="100%" x2="0%" y2="0%">
                      <stop offset="0%" stopColor="hsl(270, 80%, 50%)" />
                      <stop offset="50%" stopColor="hsl(165, 80%, 50%)" />
                      <stop offset="100%" stopColor="hsl(60, 80%, 50%)" />
                    </linearGradient>
                  </defs>

                  {/* Grid Lines — Y axis (binding, log10) */}
                  {getLogTicks(scatterBounds.yMin, scatterBounds.yMax).map((exp) => (
                    <line
                      key={`y-${exp}`}
                      x1="40"
                      y1={getScatterY(Math.pow(10, exp))}
                      x2="400"
                      y2={getScatterY(Math.pow(10, exp))}
                      className="stroke-border"
                      strokeWidth="1"
                    />
                  ))}
                  {/* Grid Lines — X axis (wtLeakage, log10) */}
                  {getLogTicks(scatterBounds.xMin, scatterBounds.xMax).map((exp) => (
                    <line
                      key={`x-${exp}`}
                      x1={getScatterX(Math.pow(10, exp))}
                      y1="20"
                      x2={getScatterX(Math.pow(10, exp))}
                      y2="170"
                      className="stroke-border"
                      strokeWidth="1"
                    />
                  ))}

                  {/* Axes */}
                  <line x1="40" y1="170" x2="400" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />
                  <line x1="40" y1="20" x2="40" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />

                  {/* X Axis Labels */}
                  {getLogTicks(scatterBounds.xMin, scatterBounds.xMax).map((exp) => (
                    <text
                      key={`xl-${exp}`}
                      x={getScatterX(Math.pow(10, exp))}
                      y="185"
                      textAnchor="middle"
                      className="text-[9px] fill-muted-foreground font-bold"
                    >
                      10<sup>{exp}</sup>
                    </text>
                  ))}
                  <text x="220" y="198" textAnchor="middle" className="text-[10px] fill-muted-foreground font-bold">
                    WT Leakage TIR (lower is better)
                  </text>

                  {/* Y Axis Labels */}
                  {getLogTicks(scatterBounds.yMin, scatterBounds.yMax).map((exp) => (
                    <text
                      key={`yl-${exp}`}
                      x="34"
                      y={getScatterY(Math.pow(10, exp)) + 3}
                      textAnchor="end"
                      className="text-[9px] fill-muted-foreground font-bold"
                    >
                      10<sup>{exp}</sup>
                    </text>
                  ))}
                  <text
                    x="12"
                    y="95"
                    textAnchor="middle"
                    transform="rotate(-90, 12, 95)"
                    className="text-[10px] fill-muted-foreground font-bold"
                  >
                    Orthogonal Binding (higher is better)
                  </text>

                  {/* Plot Scatter Points */}
                  {scatterPoints.map((pt) => {
                    const isHovered = hoveredScatterPoint?.id === pt.id;
                    return (
                      <circle
                        key={pt.id}
                        cx={getScatterX(pt.wtLeakage)}
                        cy={getScatterY(pt.binding)}
                        r={isHovered ? 6 : 2.5}
                        fill={getScatterColor(pt.access)}
                        className="stroke-card cursor-pointer transition-all duration-150"
                        strokeWidth={isHovered ? 1.5 : 0.5}
                        onMouseEnter={() =>
                          setHoveredScatterPoint({
                            ...pt,
                            x: getScatterX(pt.wtLeakage),
                            y: getScatterY(pt.binding),
                          })
                        }
                        onMouseLeave={() => setHoveredScatterPoint(null)}
                      />
                    );
                  })}

                  {/* Legend Slider */}
                  <g transform="translate(425, 20)">
                    <rect x="0" y="0" width="12" height="150" fill="url(#access-grad)" rx="2" />
                    <text x="20" y="10" className="text-[8px] fill-muted-foreground font-bold">High</text>
                    <text x="20" y="146" className="text-[8px] fill-muted-foreground font-bold">Low</text>

                    <text
                      x="36"
                      y="75"
                      textAnchor="middle"
                      transform="rotate(90, 36, 75)"
                      className="text-[9px] fill-muted-foreground font-bold"
                    >
                      Accessibility
                    </text>
                  </g>
                </svg>

                {/* Scatter Point Tooltip */}
                {hoveredScatterPoint && (
                  <div
                    className="absolute bg-neutral-900/95 text-white text-[10px] p-2.5  shadow-lg pointer-events-none z-30 flex flex-col space-y-1 border border-neutral-800"
                    style={{
                      left: `${(hoveredScatterPoint.x / 500) * 100}%`,
                      top: `${(hoveredScatterPoint.y / 200) * 100}%`,
                      transform: hoveredScatterPoint.x > 250 ? "translate(-110%, -110%)" : "translate(10%, -110%)",
                    }}
                  >
                    <div className="flex justify-between space-x-3 border-b border-neutral-800 pb-1 mb-1">
                      <span className="font-bold text-neutral-300">Candidate Details</span>
                    </div>
                    <div className="flex justify-between space-x-4">
                      <span className="text-neutral-400 font-semibold">WT Leakage TIR:</span>
                      <span className="font-bold">
                        {hoveredScatterPoint.wtLeakage.toExponential(2)}
                      </span>
                    </div>
                    <div className="flex justify-between space-x-4">
                      <span className="text-neutral-400 font-semibold">Orth. Binding TIR:</span>
                      <span className="font-bold text-red-400">
                        {hoveredScatterPoint.binding.toExponential(2)}
                      </span>
                    </div>
                    <div className="flex justify-between space-x-4">
                      <span className="text-neutral-400 font-semibold">Accessibility:</span>
                      <span className="font-bold" style={{ color: getScatterColor(hoveredScatterPoint.access) }}>
                        {hoveredScatterPoint.access.toFixed(2)}
                      </span>
                    </div>
                  </div>
                )}
              </div>
            </div>

            {/* CARD E: Top Candidates Table */}
            <div className="col-span-2 row-span-2 bg-card text-card-foreground p-6  border border-border  flex flex-col justify-between h-full min-h-[380px]">
              <div>
                <h3 className="font-bold text-card-foreground text-xs tracking-wider uppercase mb-3">
                  Top Candidates
                </h3>

                <div className="overflow-x-auto">
                  <table className="w-full text-left text-[11px] border-collapse">
                    <thead>
                      <tr className="border-b border-border text-muted-foreground font-bold">
                        <th className="py-2.5 pb-2">Rank</th>
                        <th className="py-2.5 pb-2">RBS</th>
                        <th className="py-2.5 pb-2">Spacer</th>
                        <th className="py-2.5 pb-2 text-right">Orth</th>
                        <th className="py-2.5 pb-2 text-right">WT Leak</th>
                        <th className="py-2.5 pb-2 text-right">Access</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-border/60">
                      {candidates.map((cand, idx) => {
                        const isSelected = selectedCandidateIndex === idx;
                        return (
                          <tr
                            key={idx}
                            onClick={() => setSelectedCandidateIndex(idx)}
                            className={`cursor-pointer hover:bg-muted/40 transition-colors ${isSelected
                              ? "bg-primary/5 text-primary font-bold border-l-2 border-primary"
                              : ""
                              }`}
                          >
                            <td className="py-3 font-bold">
                              {idx === 0 ? (
                                <span className="inline-flex items-center justify-center w-5 h-5 rounded-full bg-[#fef08a] dark:bg-amber-950/80 text-amber-800 dark:text-amber-300 text-[10px]">🥇</span>
                              ) : idx === 1 ? (
                                <span className="inline-flex items-center justify-center w-5 h-5 rounded-full bg-neutral-100 dark:bg-neutral-850 text-neutral-700 dark:text-neutral-300 text-[10px]">🥈</span>
                              ) : idx === 2 ? (
                                <span className="inline-flex items-center justify-center w-5 h-5 rounded-full bg-amber-100 dark:bg-amber-950/80 text-amber-800 dark:text-amber-300 text-[10px]">🥉</span>
                              ) : (
                                <span className="pl-1.5">{idx + 1}</span>
                              )}
                            </td>
                            <td className="py-3 font-mono tracking-tight text-foreground">
                              {cand.rbs}
                            </td>
                            <td className="py-3 text-muted-foreground font-medium">
                              <span className="font-mono">{cand.spacer}</span>
                              <span className="text-[9px] text-muted-foreground/80 ml-1">({cand.spacer.length}n)</span>
                            </td>
                            <td className="py-3 text-right font-bold text-primary">
                              {cand.orthScore}
                            </td>
                            <td className="py-3 text-right font-bold text-muted-foreground">
                              {cand.wtLeakage}
                            </td>
                            <td className="py-3 text-right font-bold text-foreground/80">
                              {cand.rbsAccess}
                            </td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
              </div>

              <div className="text-[10px] text-muted-foreground border-t border-border pt-3 mt-3 leading-relaxed font-semibold">
                All scores are normalized (0-1). Higher Orth Score and Fitness are better. Lower WT Leakage is better.
              </div>
            </div>

            {/* CARD F: Local Structure Inspector */}
            <div className="col-span-2 row-span-2 bg-card text-card-foreground p-6  border border-border  flex flex-col justify-between h-full min-h-[380px]">
              <div>
                <div className="flex items-center space-x-1.5 mb-4 border-b border-border pb-2">
                  <h3 className="font-bold text-card-foreground text-xs tracking-wider uppercase">
                    Structure Inspector
                  </h3>
                  <div className="relative cursor-pointer group/info">
                    <Info className="w-4 h-4 text-muted-foreground hover:text-foreground" />
                    <div className="absolute left-1/2 -translate-x-1/2 bottom-full mb-2 w-48 bg-neutral-900 text-white text-[10px] p-2.5  opacity-0 pointer-events-none group-hover/info:opacity-100 transition-opacity z-50 shadow-lg leading-relaxed">
                      Displays accessibility of RBS sequence under folding. High accessibility of RBS ensures efficient translation start.
                    </div>
                  </div>
                </div>

                <div className="flex justify-between text-[10px] font-bold text-muted-foreground mb-2">
                  <span>Sequence (5' &rarr; 3')</span>
                  <span>~45 nt window</span>
                </div>

                {/* Sequence box - Red highlight */}
                <div className="bg-muted/40 p-3.5  border border-border font-mono text-[11px] overflow-x-auto whitespace-nowrap mb-4 leading-relaxed text-muted-foreground select-all">
                  <span className="opacity-60">... GCUUU </span>
                  <span className="bg-primary/10 text-primary font-extrabold px-1 py-0.5 rounded">
                    {selectedCandidate ? selectedCandidate.rbs : ""}
                  </span>
                  <span className="font-bold text-foreground/80">
                    {" "}{selectedCandidate ? selectedCandidate.spacer : ""}ACGACAA{" "}
                  </span>
                  <span className="bg-blue-100/10 text-blue-500 font-extrabold px-1 py-0.5 rounded">
                    AUG
                  </span>
                  <span className="opacity-60"> GCUACU ...</span>
                </div>

                <div className="text-[10px] font-bold text-muted-foreground mb-2">
                  Structure (RNAfold)
                </div>

                {/* Structure box - Red highlight */}
                <div className="bg-muted/40 p-3.5  border border-border font-mono text-[11px] overflow-x-auto whitespace-nowrap leading-relaxed text-muted-foreground select-all">
                  <span className="opacity-60">... </span>
                  <span className="bg-primary/10 text-primary font-extrabold px-1 py-0.5 rounded" title="RBS Site">
                    {selectedCandidate ? (() => {
                      const flankLen = selectedCandidate.five_prime_flank ? selectedCandidate.five_prime_flank.length : 0;
                      return selectedCandidate.structure.substring(flankLen, flankLen + selectedCandidate.rbs.length);
                    })() : ""}
                  </span>
                  <span className="text-muted-foreground font-bold" title="Spacer">
                    {selectedCandidate ? (() => {
                      const flankLen = selectedCandidate.five_prime_flank ? selectedCandidate.five_prime_flank.length : 0;
                      const rbsLen = selectedCandidate.rbs.length;
                      return selectedCandidate.structure.substring(flankLen + rbsLen, flankLen + rbsLen + selectedCandidate.spacer.length);
                    })() : ""}
                  </span>
                  <span className="bg-blue-100/10 text-blue-500 font-extrabold px-1 py-0.5 rounded" title="Start Codon">
                    {selectedCandidate ? (() => {
                      const flankLen = selectedCandidate.five_prime_flank ? selectedCandidate.five_prime_flank.length : 0;
                      const rbsLen = selectedCandidate.rbs.length;
                      const spacerLen = selectedCandidate.spacer.length;
                      const start = flankLen + rbsLen + spacerLen;
                      if (selectedCandidate.five_prime_flank && start + 3 <= selectedCandidate.structure.length) {
                        return selectedCandidate.structure.substring(start, start + 3);
                      }
                      return selectedCandidate.structure.substring(selectedCandidate.structure.length - 3);
                    })() : ""}
                  </span>
                  <span className="text-muted-foreground opacity-60" title="Downstream remainder">
                    {selectedCandidate ? (() => {
                      const flankLen = selectedCandidate.five_prime_flank ? selectedCandidate.five_prime_flank.length : 0;
                      const rbsLen = selectedCandidate.rbs.length;
                      const spacerLen = selectedCandidate.spacer.length;
                      const start = flankLen + rbsLen + spacerLen + 3;
                      if (selectedCandidate.five_prime_flank && start < selectedCandidate.structure.length) {
                        return selectedCandidate.structure.substring(start);
                      }
                      return "";
                    })() : ""}
                  </span>
                  <span className="opacity-60"> .....</span>
                </div>
              </div>

              {/* Legend Footer */}
              <div className="flex items-center justify-around text-[10px] font-bold text-muted-foreground border-t border-border pt-3 mt-4">
                <div className="flex items-center space-x-1.5">
                  <span className="w-2.5 h-2.5 bg-primary rounded-sm inline-block" />
                  <span>Accessible</span>
                </div>
                <div className="flex items-center space-x-1.5">
                  <span className="w-1.5 h-1.5 bg-muted-foreground rounded-full inline-block" />
                  <span>Paired</span>
                </div>
                <div className="flex items-center space-x-1.5">
                  <span className="w-2.5 h-2.5 bg-blue-400 rounded-sm inline-block" />
                  <span>Start Codon</span>
                </div>
              </div>
            </div>
          </div>
        </div>

      </div>
    </div>
  );
}
