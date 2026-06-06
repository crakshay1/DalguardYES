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
    fitness: "0.91",
    structure: "((((((.......))),). ... ((...)) AUG .....",
  },
  {
    rbs: "GGAGGAG",
    spacer: "UACUUA",
    orthScore: "0.89",
    wtLeakage: "0.04",
    rbsAccess: "0.85",
    fitness: "0.88",
    structure: "((((((.......))))). ... ((...)) AUG .....",
  },
  {
    rbs: "AGGAGGA",
    spacer: "ACGUUA",
    orthScore: "0.86",
    wtLeakage: "0.05",
    rbsAccess: "0.82",
    fitness: "0.84",
    structure: "(((((.........))))). ... ((...)) AUG .....",
  },
  {
    rbs: "GAGAGGA",
    spacer: "UACGCU",
    orthScore: "0.83",
    wtLeakage: "0.06",
    rbsAccess: "0.78",
    fitness: "0.80",
    structure: "(((((.........))))). ... ((...)) AUG .....",
  },
  {
    rbs: "GGAAAGA",
    spacer: "ACGUUU",
    orthScore: "0.80",
    wtLeakage: "0.07",
    rbsAccess: "0.76",
    fitness: "0.77",
    structure: "((((...........)))). ... ((...)) AUG .....",
  },
];

const INITIAL_FITNESS_DATA = (() => {
  const data = [];
  for (let i = 0; i <= 100; i++) {
    const bestNoise = Math.sin(i / 6) * 0.006 + (i % 7 === 0 ? 0.004 : -0.004);
    const avgNoise = Math.cos(i / 10) * 0.005 + (i % 5 === 0 ? 0.002 : -0.002);

    const bestVal = Math.min(1.0, 0.18 + 0.60 * Math.pow(i / 100, 0.35) + bestNoise);
    const avgVal = Math.min(bestVal, 0.10 + 0.48 * Math.pow(i / 100, 0.5) + avgNoise);

    data.push({
      generation: i,
      best: parseFloat(bestVal.toFixed(3)),
      avg: parseFloat(avgVal.toFixed(3)),
    });
  }
  return data;
})();

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
  fitness: string;
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
  const [sequence, setSequence] = useState("");
  const [fileName, setFileName] = useState("");

  // Input states
  const [orthogonalAntiSD, setOrthogonalAntiSD] = useState("");
  const [wtAntiSD, setWtAntiSD] = useState("");
  const [cdsStart, setCdsStart] = useState("");
  const [targetExpression, setTargetExpression] = useState("High");

  // Selection state
  const [selectedCandidateIndex, setSelectedCandidateIndex] = useState(0);

  // Copy success indicator
  const [copiedField, setCopiedField] = useState<string | null>(null);

  // Optimization run simulation
  const [isOptimizing, setIsOptimizing] = useState(false);
  const [optProgress, setOptProgress] = useState(0);

  // Stateful datasets
  const [candidates, setCandidates] = useState<Candidate[]>(INITIAL_CANDIDATES);
  const [fitnessData, setFitnessData] = useState<FitnessPoint[]>(INITIAL_FITNESS_DATA);
  const [scatterPoints, setScatterPoints] = useState<ScatterPoint[]>(INITIAL_SCATTER_POINTS);

  // Parse uploaded file sequence
  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    setFileName(file.name);
    const reader = new FileReader();
    reader.onload = (event) => {
      const text = event.target?.result as string;
      if (!text) return;
      const parsed = parseFASTA(text, file.name);
      setSequence(parsed.sequence);
    };
    reader.readAsText(file);
  };

  // Populate inputs when sequence file is uploaded
  useEffect(() => {
    if (sequence) {
      if (sequence.length >= 9) {
        setOrthogonalAntiSD(sequence.substring(0, 9));
      }
      if (sequence.length >= 18) {
        setWtAntiSD(sequence.substring(9, 18));
      }
      if (sequence.length >= 27) {
        setCdsStart(sequence.substring(18, 27) + "...");
      }
    }
  }, [sequence]);

  const handleCopy = (text: string, fieldId: string) => {
    navigator.clipboard.writeText(text);
    setCopiedField(fieldId);
    setTimeout(() => setCopiedField(null), 1500);
  };

  const startOptimization = () => {
    if (isOptimizing) return;
    setIsOptimizing(true);
    setOptProgress(0);

    const interval = setInterval(() => {
      setOptProgress((prev) => {
        if (prev >= 100) {
          clearInterval(interval);
          setIsOptimizing(false);
          // Shuffle candidate ranks for simulation feedback
          setSelectedCandidateIndex((prevIndex) => (prevIndex + 1) % candidates.length);
          return 100;
        }
        return prev + 10;
      });
    }, 150);
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
            if (parsed.inputs.orthogonalAntiSD) setOrthogonalAntiSD(parsed.inputs.orthogonalAntiSD);
            if (parsed.inputs.wtAntiSD) setWtAntiSD(parsed.inputs.wtAntiSD);
            if (parsed.inputs.cdsStart) setCdsStart(parsed.inputs.cdsStart);
            if (parsed.inputs.targetExpression) setTargetExpression(parsed.inputs.targetExpression);
          }

          if (Array.isArray(parsed.candidates)) {
            setCandidates(parsed.candidates);
            setSelectedCandidateIndex(0);
          }

          if (Array.isArray(parsed.fitnessData)) {
            setFitnessData(parsed.fitnessData);
          }

          if (Array.isArray(parsed.scatterPoints)) {
            setScatterPoints(parsed.scatterPoints);
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
        } else if (result.type === "fitness") {
          setFitnessData(result.data);
          alert(`Loaded ${result.data.length} generations of fitness data from CSV!`);
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
      const upstream = "GCTTT";
      const rbsSeq = selectedCandidate.rbs.replace(/U/g, "T");
      const spacerSeq = selectedCandidate.spacer.replace(/U/g, "T");
      const startCodon = "ATG";

      // Calculate remaining CDS sequence from cdsStart field
      const cdsRemaining = cdsStart
        .replace(/\./g, "")
        .replace(/U/g, "T")
        .toUpperCase();

      const cdsPart = cdsRemaining.startsWith("ATG") ? cdsRemaining.substring(3) : cdsRemaining;
      const computedSeq = upstream + rbsSeq + spacerSeq + startCodon + cdsPart;

      const newAnnotations = [
        { name: "Upstream", start: 0, end: upstream.length, direction: 1, color: "#6b7280" },
        { name: "RBS Site", start: upstream.length, end: upstream.length + rbsSeq.length, direction: 1, color: "#dc2626" },
        { name: "Spacer", start: upstream.length + rbsSeq.length, end: upstream.length + rbsSeq.length + spacerSeq.length, direction: 1, color: "#eab308" },
        { name: "Start Codon", start: upstream.length + rbsSeq.length + spacerSeq.length, end: upstream.length + rbsSeq.length + spacerSeq.length + 3, direction: 1, color: "#3b82f6" }
      ];

      if (cdsPart.length > 0) {
        newAnnotations.push({
          name: "Coding Sequence (CDS)",
          start: upstream.length + rbsSeq.length + spacerSeq.length + 3,
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
  }, [selectedCandidate, cdsStart]);

  const getFitnessX = (gen: number) => 40 + (gen / 100) * 440;
  const getFitnessY = (fit: number) => 170 - fit * 140;

  const bestD = useMemo(() => {
    return "M " + fitnessData.map((d) => `${getFitnessX(d.generation)},${getFitnessY(d.best)}`).join(" L ");
  }, [fitnessData]);

  const avgD = useMemo(() => {
    return "M " + fitnessData.map((d) => `${getFitnessX(d.generation)},${getFitnessY(d.avg)}`).join(" L ");
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
    const rect = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const plotWidth = 440;
    const startX = 40;
    const pct = Math.max(0, Math.min(1, (x - startX) / plotWidth));
    const gen = Math.round(pct * 100);

    const dataPoint = fitnessData[gen];
    if (dataPoint) {
      setHoveredFitnessData({
        ...dataPoint,
        x: getFitnessX(gen),
        y: getFitnessY(dataPoint.best),
      });
    }
  };

  const handleFitnessMouseLeave = () => {
    setHoveredFitnessData(null);
  };

  const getScatterX = (leakage: number) => {
    const logVal = Math.log10(leakage);
    const minLog = -4.0;
    const maxLog = 0.0;
    const pct = (logVal - minLog) / (maxLog - minLog);
    return 40 + pct * 360;
  };

  const getScatterY = (binding: number) => {
    return 170 - binding * 150;
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

            {/* CARD A: Inputs */}
            <div className="col-span-3 bg-card text-card-foreground p-6  border border-border  flex flex-col justify-between">
              <div className="flex items-center justify-between mb-4 border-b border-border pb-3">
                <h2 className="font-bold text-card-foreground/80 text-xs tracking-wider uppercase">
                  Designer Inputs
                </h2>

                {/* File Upload Trigger */}
                <label className="flex items-center space-x-1.5 text-xs text-primary hover:text-primary/80 font-semibold cursor-pointer transition-colors bg-primary/5 px-2.5 py-1 rounded-md border border-primary/10">
                  <Upload className="w-3.5 h-3.5" />
                  <span>{fileName ? fileName : "Upload sequence"}</span>
                  <input
                    type="file"
                    accept=".fasta,.txt,.seq"
                    className="hidden"
                    onChange={handleFileUpload}
                  />
                </label>
              </div>

              <div className="grid grid-cols-2 gap-4">
                {/* Input 1 */}
                <div className="space-y-1">
                  <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                    Orthogonal anti-SD
                  </label>
                  <div className="relative flex items-center">
                    <input
                      type="text"
                      value={orthogonalAntiSD}
                      onChange={(e) => setOrthogonalAntiSD(e.target.value.toUpperCase())}
                      className="w-full bg-muted/40 border border-border  px-3 py-2 text-sm font-mono focus:outline-none focus:ring-1 focus:ring-primary pr-10 uppercase text-foreground"
                    />
                    <button
                      onClick={() => handleCopy(orthogonalAntiSD, "orth")}
                      className="absolute right-2 text-muted-foreground hover:text-foreground p-1"
                    >
                      {copiedField === "orth" ? (
                        <Check className="w-4 h-4 text-primary animate-in fade-in zoom-in-50 duration-200" />
                      ) : (
                        <Copy className="w-4 h-4" />
                      )}
                    </button>
                  </div>
                </div>

                {/* Input 2 */}
                <div className="space-y-1">
                  <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                    WT anti-SD
                  </label>
                  <div className="relative flex items-center">
                    <input
                      type="text"
                      value={wtAntiSD}
                      onChange={(e) => setWtAntiSD(e.target.value.toUpperCase())}
                      className="w-full bg-muted/40 border border-border  px-3 py-2 text-sm font-mono pr-10 focus:outline-none focus:ring-1 focus:ring-primary uppercase text-foreground"
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

                {/* Input 3 */}
                <div className="space-y-1">
                  <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                    CDS start
                  </label>
                  <div className="relative flex items-center">
                    <input
                      type="text"
                      value={cdsStart}
                      onChange={(e) => setCdsStart(e.target.value)}
                      className="w-full bg-muted/40 border border-border  px-3 py-2 text-sm font-mono pr-10 focus:outline-none focus:ring-1 focus:ring-primary text-foreground"
                    />
                    <button
                      onClick={() => handleCopy(cdsStart, "cds")}
                      className="absolute right-2 text-muted-foreground hover:text-foreground p-1"
                    >
                      {copiedField === "cds" ? (
                        <Check className="w-4 h-4 text-primary animate-in fade-in zoom-in-50 duration-200" />
                      ) : (
                        <Copy className="w-4 h-4" />
                      )}
                    </button>
                  </div>
                </div>

                {/* Input 4 */}
                <div className="space-y-1">
                  <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                    Target expression
                  </label>
                  <div className="relative">
                    <select
                      value={targetExpression}
                      onChange={(e) => setTargetExpression(e.target.value)}
                      className="w-full bg-muted/40 border border-border  px-3 py-2 text-sm appearance-none focus:outline-none focus:ring-1 focus:ring-primary cursor-pointer text-foreground"
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
                    Tracks the maximum (Best) and mean (Average) fitness scores of RBS designs over 100 generations of the Genetic Algorithm.
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
                  {[0, 20, 40, 60, 80, 100].map((gen) => (
                    <line
                      key={gen}
                      x1={getFitnessX(gen)}
                      y1="20"
                      x2={getFitnessX(gen)}
                      y2="170"
                      className="stroke-border"
                      strokeWidth="1"
                    />
                  ))}

                  {/* Axes */}
                  <line x1="40" y1="170" x2="480" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />
                  <line x1="40" y1="20" x2="40" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />

                  {/* Axis Labels */}
                  {[0, 20, 40, 60, 80, 100].map((gen) => (
                    <text
                      key={gen}
                      x={getFitnessX(gen)}
                      y="185"
                      textAnchor="middle"
                      className="text-[9px] fill-muted-foreground font-bold"
                    >
                      {gen}
                    </text>
                  ))}
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

                  {/* Grid Lines */}
                  {[0, 0.2, 0.4, 0.6, 0.8, 1.0].map((val) => (
                    <line
                      key={val}
                      x1="40"
                      y1={getScatterY(val)}
                      x2="400"
                      y2={getScatterY(val)}
                      className="stroke-border"
                      strokeWidth="1"
                    />
                  ))}
                  {[-4, -3, -2, -1, 0].map((exp) => {
                    const wtLeak = Math.pow(10, exp);
                    return (
                      <line
                        key={exp}
                        x1={getScatterX(wtLeak)}
                        y1="20"
                        x2={getScatterX(wtLeak)}
                        y2="170"
                        className="stroke-border"
                        strokeWidth="1"
                      />
                    );
                  })}


                  {/* Axes */}
                  <line x1="40" y1="170" x2="400" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />
                  <line x1="40" y1="20" x2="40" y2="170" className="stroke-muted-foreground/50" strokeWidth="1.5" />

                  {/* Labels */}
                  {[-4, -3, -2, -1, 0].map((exp) => {
                    const wtLeak = Math.pow(10, exp);
                    return (
                      <text
                        key={exp}
                        x={getScatterX(wtLeak)}
                        y="185"
                        textAnchor="middle"
                        className="text-[9px] fill-muted-foreground font-bold"
                      >
                        10<sup>{exp}</sup>
                      </text>
                    );
                  })}
                  <text x="220" y="198" textAnchor="middle" className="text-[10px] fill-muted-foreground font-bold">
                    WT Leakage (lower is better)
                  </text>

                  {[0, 0.2, 0.4, 0.6, 0.8, 1.0].map((val) => (
                    <text
                      key={val}
                      x="30"
                      y={getScatterY(val) + 3}
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
                      <span className="text-neutral-400 font-semibold">WT Leakage:</span>
                      <span className="font-bold">{hoveredScatterPoint.wtLeakage.toExponential(2)}</span>
                    </div>
                    <div className="flex justify-between space-x-4">
                      <span className="text-neutral-400 font-semibold">Orthogonal Binding:</span>
                      <span className="font-bold text-red-400">{hoveredScatterPoint.binding.toFixed(2)}</span>
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
                        <th className="py-2.5 pb-2 text-right">Fitness</th>
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
                            <td className="py-3 text-right font-extrabold text-foreground">
                              {cand.fitness}
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
                  <span className="bg-primary/10 text-primary font-extrabold px-1 py-0.5 rounded">
                    {selectedCandidate ? selectedCandidate.structure.substring(0, selectedCandidate.rbs.length) : ""}
                  </span>
                  <span className="text-muted-foreground font-bold">
                    {selectedCandidate ? selectedCandidate.structure.substring(
                      selectedCandidate.rbs.length,
                      selectedCandidate.structure.length - 3
                    ) : ""}
                  </span>
                  <span className="bg-blue-100/10 text-blue-500 font-extrabold px-1 py-0.5 rounded">
                    {selectedCandidate ? selectedCandidate.structure.substring(selectedCandidate.structure.length - 3) : ""}
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

            {/* CARD G: Design Insights */}
            <div className="col-span-2 bg-card text-card-foreground p-6  border border-border  flex flex-col justify-between h-[180px]">
              <div className="flex items-center space-x-1.5 mb-2">
                <h3 className="font-bold text-card-foreground text-xs tracking-wider uppercase">
                  Design Insights
                </h3>
                <div className="relative cursor-pointer group/info">
                  <Info className="w-4 h-4 text-muted-foreground hover:text-foreground" />
                  <div className="absolute left-1/2 -translate-x-1/2 bottom-full mb-2 w-48 bg-neutral-900 text-white text-[10px] p-2.5  opacity-0 pointer-events-none group-hover/info:opacity-100 transition-opacity z-50 shadow-lg leading-relaxed">
                    Shows relative contribution of each parameter to final fitness scores calculated by RiboGuard AI model.
                  </div>
                </div>
              </div>

              <div className="space-y-2.5">
                {[
                  { label: "RBS Accessibility", value: 0.46 },
                  { label: "Orthogonal Binding", value: 0.31 },
                  { label: "Spacer Length", value: 0.15 },
                  { label: "WT Leakage", value: 0.08 },
                ].map((feat) => (
                  <div key={feat.label} className="space-y-0.5">
                    <div className="flex justify-between text-[11px]">
                      <span className="font-bold text-muted-foreground">{feat.label}</span>
                      <span className="font-extrabold text-foreground">{feat.value.toFixed(2)}</span>
                    </div>
                    <div className="h-2 w-full bg-muted border border-border/30 rounded-full overflow-hidden">
                      <div
                        className="h-full bg-primary rounded-full transition-all duration-500"
                        style={{ width: `${feat.value * 100}%` }}
                      />
                    </div>
                  </div>
                ))}
              </div>
            </div>

            {/* CARD H: Validation Status */}
            <div className="col-span-2 bg-card text-card-foreground p-6  border border-border  flex flex-col justify-between h-[180px]">
              <div className="flex items-center space-x-1.5 mb-2 border-b border-border pb-2">
                <h3 className="font-bold text-card-foreground text-xs tracking-wider uppercase">
                  Validation Status
                </h3>
                <div className="relative cursor-pointer group/info">
                  <Info className="w-4 h-4 text-muted-foreground hover:text-foreground" />
                  <div className="absolute left-1/2 -translate-x-1/2 bottom-full mb-2 w-48 bg-neutral-900 text-white text-[10px] p-2.5  opacity-0 pointer-events-none group-hover/info:opacity-100 transition-opacity z-50 shadow-lg leading-relaxed">
                    Secondary folding checkers using standard toolkits to double-check design reliability.
                  </div>
                </div>
              </div>

              <div className="space-y-2.5">
                {[
                  { label: "RNAfold", value: "pass" },
                  { label: "RNAduplex", value: "strong orthogonal match" },
                  { label: "IntaRNA", value: "low WT interaction" },
                ].map((val) => (
                  <div key={val.label} className="flex items-center space-x-2.5">
                    <div className="w-5 h-5 rounded-full bg-primary/5 flex items-center justify-center border border-primary/10">
                      <CheckCircle2 className="w-3.5 h-3.5 text-primary" />
                    </div>
                    <div className="text-xs">
                      <span className="font-bold text-muted-foreground mr-1.5">{val.label}:</span>
                      <span className="font-extrabold text-primary">{val.value}</span>
                    </div>
                  </div>
                ))}
              </div>
            </div>

          </div>
        </div>

      </div>
    </div>
  );
}
