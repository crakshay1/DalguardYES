import { useState, useEffect } from "react";
import { ThemeToggle } from "../components/ThemeToggle";
import { parseFASTA } from "../utils/dataIngestion";
import type { MetaFunction } from "react-router";
import { useNavigate } from "react-router";
import {
    Layers,
    Search,
    Upload,
    Plus,
    Trash2,
    Activity,
    Calculator,
    Dna,
    Scissors,
    FileText,
} from "lucide-react";

export const meta: MetaFunction = () => {
    return [{ title: "RiboGuard AI - Structure Inspector" }, { name: "Inspector", content: "Inspect plasmid maps and linear sequences" }];
};

interface Annotation {
    name: string;
    start: number;
    end: number;
    direction: number; // 1 for forward, -1 for reverse
    color: string;
}

// Helper to calculate GC content
const calculateGCContent = (seq: string) => {
    if (!seq) return 0;
    const gcCount = (seq.match(/[GC]/gi) || []).length;
    return (gcCount / seq.length) * 100;
};

// Helper to calculate molecular weight (g/mol)
const calculateMolecularWeight = (seq: string) => {
    if (!seq) return 0;
    const isRNA = seq.includes("U");
    const a = (seq.match(/A/gi) || []).length;
    const c = (seq.match(/C/gi) || []).length;
    const g = (seq.match(/G/gi) || []).length;
    const t = (seq.match(/[TU]/gi) || []).length;

    if (isRNA) {
        return a * 329.2 + t * 306.2 + c * 305.2 + g * 345.2 + 159.0;
    } else {
        return a * 313.21 + t * 304.2 + c * 289.18 + g * 329.21 + 79.0;
    }
};



// Codon translation table
const CODON_TABLE: Record<string, string> = {
    ATT: "I", ATC: "I", ATA: "I",
    CTT: "L", CTC: "L", CTA: "L", CTG: "L", TTA: "L", TTG: "L",
    GTT: "V", GTC: "V", GTA: "V", GTG: "V",
    TTT: "F", TTC: "F",
    ATG: "M",
    TGT: "C", TGC: "C",
    GCT: "A", GCC: "A", GCA: "A", GCG: "A",
    GGT: "G", GGC: "G", GGA: "G", GGG: "G",
    CCT: "P", CCC: "P", CCA: "P", CCG: "P",
    ACT: "T", ACC: "T", ACA: "T", ACG: "T",
    TCT: "S", TCC: "S", TCA: "S", TCG: "S", AGT: "S", AGC: "S",
    TAT: "Y", TAC: "Y",
    TGG: "W",
    CAA: "Q", CAG: "Q",
    AAT: "N", AAC: "N",
    CAT: "H", CAC: "H",
    GAA: "E", GAG: "E",
    GAT: "D", GAC: "D",
    AAA: "K", AAG: "K",
    CGT: "R", CGC: "R", CGA: "R", CGG: "R", AGA: "R", AGG: "R",
    TAA: "*", TAG: "*", TGA: "*"
};

// Translate sequence
const translateSequence = (seq: string) => {
    if (!seq) return "";
    const dna = seq.replace(/U/g, "T").toUpperCase();

    // Find the first start codon 'ATG'
    const startIdx = dna.indexOf("ATG");
    if (startIdx === -1) {
        // Translate from index 0 if no start codon found
        let protein = "";
        for (let i = 0; i <= dna.length - 3; i += 3) {
            const codon = dna.substring(i, i + 3);
            protein += CODON_TABLE[codon] || "?";
        }
        return protein;
    }

    let protein = "";
    for (let i = startIdx; i <= dna.length - 3; i += 3) {
        const codon = dna.substring(i, i + 3);
        const aa = CODON_TABLE[codon] || "?";
        protein += aa;
        if (aa === "*") break; // Stop codon
    }
    return protein;
};

// Standard restriction enzymes list
const RESTRICTION_ENZYMES = [
    { name: "EcoRI", site: "GAATTC" },
    { name: "BamHI", site: "GGATCC" },
    { name: "HindIII", site: "AAGCTT" },
    { name: "XhoI", site: "CTCGAG" },
    { name: "NcoI", site: "CCATGG" },
    { name: "SacI", site: "GAGCTC" },
    { name: "KpnI", site: "GGTACC" },
    { name: "SmaI", site: "CCCGGG" }
];

// Find restriction cut sites
const findRestrictionSites = (seq: string) => {
    if (!seq) return [];
    const dna = seq.replace(/U/g, "T").toUpperCase();
    const results: { name: string; site: string; positions: number[] }[] = [];

    RESTRICTION_ENZYMES.forEach((enzyme) => {
        const positions: number[] = [];
        let idx = dna.indexOf(enzyme.site);
        while (idx !== -1) {
            positions.push(idx + 1); // 1-indexed cut position
            idx = dna.indexOf(enzyme.site, idx + 1);
        }
        if (positions.length > 0) {
            results.push({ name: enzyme.name, site: enzyme.site, positions });
        }
    });
    return results;
};

export default function Inspector() {
    const navigate = useNavigate();
    const [seqName, setSeqName] = useState("RiboGuard_Design_01");
    const [SeqVizComponent, setSeqVizComponent] = useState<any>(null);
    const [sequence, setSequence] = useState(
        "ATGCTAGCGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    );
    const [structure, setStructure] = useState(
        ""
    );
    const [viewerType, setViewerType] = useState<"linear" | "circular" | "both">("linear");
    const [hoveredBaseIndex, setHoveredBaseIndex] = useState<number | null>(null);

    // Custom annotations state
    const [annotations, setAnnotations] = useState<Annotation[]>([
        { name: "Start Codon", start: 0, end: 3, direction: 1, color: "#3b82f6" },
        { name: "RBS Site", start: 12, end: 20, direction: 1, color: "#ef4444" },
        { name: "Spacer", start: 20, end: 35, direction: 1, color: "#eab308" },
    ]);

    // Load dynamic client-side SeqViz and restore computed sequence from localStorage
    useEffect(() => {
        import("seqviz")
            .then((mod) => {
                setSeqVizComponent(() => mod.SeqViz);
            })
            .catch((err) => {
                console.error("Failed to load seqviz library:", err);
            });

        // Load synced sequence and structure from design step
        const stored = localStorage.getItem("riboguard_computed_sequence");
        if (stored) {
            try {
                const parsed = JSON.parse(stored);
                if (parsed.sequence) {
                    setSequence(parsed.sequence);
                    setSeqName(parsed.name || "Computed Sequence");
                    if (parsed.structure) {
                        setStructure(parsed.structure);
                    }
                    if (parsed.annotations) {
                        setAnnotations(parsed.annotations);
                    }
                }
            } catch (e) {
                console.error("Error parsing computed sequence from localStorage", e);
            }
        }
    }, []);

    // New annotation form inputs
    const [newAnnName, setNewAnnName] = useState("");
    const [newAnnStart, setNewAnnStart] = useState(0);
    const [newAnnEnd, setNewAnnEnd] = useState(10);
    const [newAnnDirection, setNewAnnDirection] = useState(1);
    const [newAnnColor, setNewAnnColor] = useState("#ef4444");

    // File parse upload
    const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;

        const reader = new FileReader();
        reader.onload = (event) => {
            const text = event.target?.result as string;
            if (!text) return;
            const parsed = parseFASTA(text, file.name);
            setSeqName(parsed.name);
            setSequence(parsed.sequence);
        };
        reader.readAsText(file);
    };

    const addAnnotation = () => {
        if (!newAnnName.trim()) return;
        setAnnotations([
            ...annotations,
            {
                name: newAnnName,
                start: Number(newAnnStart),
                end: Number(newAnnEnd),
                direction: newAnnDirection,
                color: newAnnColor,
            },
        ]);
        // reset form
        setNewAnnName("");
    };

    const removeAnnotation = (index: number) => {
        setAnnotations(annotations.filter((_, i) => i !== index));
    };

    // Derived properties metrics
    const totalLen = sequence.length;
    const gcContent = calculateGCContent(sequence);
    const molWeight = calculateMolecularWeight(sequence);
    const restrictionSites = findRestrictionSites(sequence);
    const translation = translateSequence(sequence);

    // Get clean, aligned structure that matches the sequence length
    const getAlignedStructure = () => {
        const cleanStr = structure.replace(/[^()]/g, ".");
        if (cleanStr.length === sequence.length) return cleanStr;
        if (cleanStr.length < sequence.length) {
            return cleanStr + ".".repeat(sequence.length - cleanStr.length);
        }
        return cleanStr.substring(0, sequence.length);
    };
    const alignedStructure = getAlignedStructure();

    // Map base-pairs from structure
    const getBasePairs = (db: string) => {
        const pairs: Record<number, number> = {};
        const stack: number[] = [];
        for (let i = 0; i < db.length; i++) {
            if (db[i] === "(") {
                stack.push(i);
            } else if (db[i] === ")") {
                if (stack.length > 0) {
                    const openIdx = stack.pop()!;
                    pairs[openIdx] = i;
                    pairs[i] = openIdx;
                }
            }
        }
        return pairs;
    };
    const basePairs = getBasePairs(alignedStructure);

    // Get nucleotide frequencies
    const getFrequencies = () => {
        if (!sequence) return { A: 0, T: 0, G: 0, C: 0, A_pct: 0, T_pct: 0, G_pct: 0, C_pct: 0 };
        const a = (sequence.match(/A/gi) || []).length;
        const c = (sequence.match(/C/gi) || []).length;
        const g = (sequence.match(/G/gi) || []).length;
        const t = (sequence.match(/[TU]/gi) || []).length;
        return {
            A: a,
            T: t,
            G: g,
            C: c,
            A_pct: totalLen > 0 ? (a / totalLen) * 100 : 0,
            T_pct: totalLen > 0 ? (t / totalLen) * 100 : 0,
            G_pct: totalLen > 0 ? (g / totalLen) * 100 : 0,
            C_pct: totalLen > 0 ? (c / totalLen) * 100 : 0
        };
    };
    const freqs = getFrequencies();

    // Find annotation for a specific index position
    const getAnnotationAt = (idx: number) => {
        return annotations.find(ann => idx >= ann.start && idx < ann.end);
    };

    // Split sequence indices into horizontal wrapped chunks of 40 bases
    const chunkLength = 40;
    const chunkCount = Math.ceil(sequence.length / chunkLength);
    const chunks: number[][] = [];
    for (let c = 0; c < chunkCount; c++) {
        const chunk: number[] = [];
        for (let i = c * chunkLength; i < Math.min((c + 1) * chunkLength, sequence.length); i++) {
            chunk.push(i);
        }
        chunks.push(chunk);
    }

    return (
        <div className="bg-background text-foreground min-h-screen flex font-sans overflow-hidden">

            {/* Sidebar - Inherits from global theme */}
            <aside className="w-64 bg-sidebar text-sidebar-foreground flex flex-col justify-between shrink-0 border-r border-sidebar-border select-none h-screen">
                <div>
                    {/* Logo Header */}
                    <div className="p-6 flex items-center space-x-3 border-b border-sidebar-border/60">
                        <div className="bg-sidebar-primary p-1.5 rounded-lg text-sidebar-primary-foreground">
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
                            { name: "Design", icon: Layers, path: "/" },
                            { name: "Structure Inspector", icon: Search, path: "/inspector", active: true },
                        ].map((item) => (
                            <button
                                key={item.name}
                                onClick={() => !item.disabled && item.path && navigate(item.path)}
                                className="w-full flex items-center space-x-3 px-4 py-3 rounded-lg text-sm font-semibold transition-all duration-150 text-sidebar-foreground/70 hover:bg-sidebar-accent hover:text-sidebar-accent-foreground"
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
                            Structure Inspector
                        </h1>
                        <p className="text-xs text-muted-foreground mt-1 font-medium">
                            Analyze custom DNA/RNA/Protein sequences using LatchBio SeqViz visualization.
                        </p>
                    </div>

                </header>

                {/* Workspace Grid */}
                <div className="flex-1 overflow-y-auto p-8 space-y-6 bg-muted/20">
                    <div className="grid grid-cols-12 gap-6 max-w-[1400px] mx-auto items-stretch">

                        {/* Left Controls Card (4 columns) */}
                        <div className="col-span-4 bg-card text-card-foreground p-6 border border-border flex flex-col justify-between space-y-4">
                            <div>
                                <div className="flex items-center justify-between border-b border-border pb-3 mb-4">
                                    <h2 className="font-bold text-card-foreground/80 text-xs tracking-wider uppercase">
                                        Sequence Input & Options
                                    </h2>
                                    <label className="flex items-center space-x-1.5 text-xs text-primary hover:text-primary/80 font-semibold cursor-pointer transition-colors bg-primary/5 px-2.5 py-1 rounded-md border border-primary/10">
                                        <Upload className="w-3.5 h-3.5" />
                                        <span>Upload FASTA</span>
                                        <input
                                            type="file"
                                            accept=".fasta,.txt,.fa,.seq,.gb"
                                            className="hidden"
                                            onChange={handleFileUpload}
                                        />
                                    </label>
                                </div>

                                <div className="space-y-4">
                                    {/* Sequence Name */}
                                    <div className="space-y-1">
                                        <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                                            Sequence Name
                                        </label>
                                        <input
                                            type="text"
                                            value={seqName}
                                            onChange={(e) => setSeqName(e.target.value)}
                                            className="w-full bg-muted/40 border border-border px-3 py-2 text-sm font-semibold text-foreground focus:outline-none focus:ring-1 focus:ring-primary"
                                        />
                                    </div>

                                    {/* Raw Sequence input */}
                                    <div className="space-y-1">
                                        <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                                            Biological Sequence (FASTA/DNA/RNA/Protein)
                                        </label>
                                        <textarea
                                            value={sequence}
                                            onChange={(e) => setSequence(e.target.value.toUpperCase())}
                                            rows={5}
                                            className="w-full bg-muted/40 border border-border px-3 py-2 text-xs font-mono uppercase text-foreground focus:outline-none focus:ring-1 focus:ring-primary leading-normal resize-none"
                                        />
                                        <span className="text-[10px] text-muted-foreground font-semibold block text-right">
                                            Length: {sequence.length} bp
                                        </span>
                                    </div>

                                    {/* Dot-Bracket Structure input */}
                                    <div className="space-y-1">
                                        <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                                            Secondary Structure (Dot-Bracket)
                                        </label>
                                        <textarea
                                            value={structure}
                                            onChange={(e) => setStructure(e.target.value)}
                                            rows={3}
                                            className="w-full bg-muted/40 border border-border px-3 py-2 text-xs font-mono text-foreground focus:outline-none focus:ring-1 focus:ring-primary leading-normal resize-none"
                                            placeholder=""
                                        />
                                        <span className="text-[10px] text-muted-foreground font-semibold block text-right">
                                            Aligned Length: {alignedStructure.length} bp
                                        </span>
                                    </div>

                                    {/* Viewer Toggle options */}
                                    <div className="space-y-1">
                                        <label className="block text-[10px] font-bold text-muted-foreground uppercase tracking-wide">
                                            Viewer Presentation
                                        </label>
                                        <div className="grid grid-cols-3 gap-2">
                                            {["linear", "circular", "both"].map((opt) => (
                                                <button
                                                    key={opt}
                                                    onClick={() => setViewerType(opt as any)}
                                                    className="py-1.5 px-3 rounded-lg text-xs font-bold capitalize transition-all border bg-card border-border hover:bg-muted/30"
                                                >
                                                    {opt}
                                                </button>
                                            ))}
                                        </div>
                                    </div>
                                </div>
                            </div>

                            {/* Dynamic Annotations Manager */}
                            <div className="border-t border-border pt-4">
                                <h3 className="font-bold text-[10px] uppercase text-muted-foreground tracking-wide mb-2">
                                    Features & Annotations
                                </h3>
                                <div className="max-h-[120px] overflow-y-auto space-y-1.5 mb-3 border border-border/40 p-2 rounded-lg bg-muted/10">
                                    {annotations.length === 0 ? (
                                        <div className="text-[11px] text-muted-foreground text-center py-4 font-medium">No annotations added yet</div>
                                    ) : (
                                        annotations.map((ann, idx) => (
                                            <div key={idx} className="flex items-center justify-between bg-muted/40 p-1.5 px-2 rounded-md border border-border/40 text-xs">
                                                <div className="flex items-center space-x-2">
                                                    <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: ann.color }} />
                                                    <span className="font-bold text-foreground">{ann.name}</span>
                                                    <span className="text-[10px] text-muted-foreground">({ann.start}-{ann.end} bp)</span>
                                                </div>
                                                <button onClick={() => removeAnnotation(idx)} className="text-muted-foreground hover:text-destructive">
                                                    <Trash2 className="w-3.5 h-3.5" />
                                                </button>
                                            </div>
                                        ))
                                    )}
                                </div>

                                {/* Add Annotation Form */}
                                <div className="grid grid-cols-12 gap-2 text-xs">
                                    <input
                                        type="text"
                                        placeholder="Name"
                                        value={newAnnName}
                                        onChange={(e) => setNewAnnName(e.target.value)}
                                        className="col-span-5 bg-muted/40 border border-border rounded-md px-2 py-1 focus:outline-none"
                                    />
                                    <input
                                        type="number"
                                        placeholder="Start"
                                        value={newAnnStart}
                                        onChange={(e) => setNewAnnStart(Number(e.target.value))}
                                        className="col-span-2 bg-muted/40 border border-border rounded-md px-1 py-1 focus:outline-none"
                                    />
                                    <input
                                        type="number"
                                        placeholder="End"
                                        value={newAnnEnd}
                                        onChange={(e) => setNewAnnEnd(Number(e.target.value))}
                                        className="col-span-2 bg-muted/40 border border-border rounded-md px-1 py-1 focus:outline-none"
                                    />
                                    <input
                                        type="color"
                                        value={newAnnColor}
                                        onChange={(e) => setNewAnnColor(e.target.value)}
                                        className="col-span-2 h-7 w-full border border-border rounded-md cursor-pointer bg-transparent"
                                    />
                                    <button
                                        onClick={addAnnotation}
                                        className="col-span-1 flex items-center justify-center bg-primary text-primary-foreground rounded-md hover:bg-primary/95"
                                    >
                                        <Plus className="w-4 h-4" />
                                    </button>
                                </div>
                            </div>
                        </div>

                        {/* Right Content Box (8 columns: Visualizer + Secondary Structure Correlation + Analytics Stack) */}
                        <div className="col-span-8 flex flex-col space-y-6">

                            {/* SeqViz Viewer Card */}
                            <div className="bg-card text-card-foreground p-6 border border-border flex flex-col justify-between items-stretch min-h-[550px]">
                                <div className="flex items-center justify-between border-b border-border pb-3 mb-4 shrink-0">
                                    <h2 className="font-bold text-card-foreground/80 text-xs tracking-wider uppercase">
                                        Sequence Visualization
                                    </h2>
                                    <span className="text-[10px] text-muted-foreground font-semibold px-2 py-0.5 bg-muted/60 border border-border">
                                        SeqViz Renderer
                                    </span>
                                </div>

                                <div className="flex-1 bg-white dark:bg-[#1a2024]/40 overflow-hidden border border-border p-4 relative min-h-[480px]">
                                    {SeqVizComponent && sequence ? (
                                        <SeqVizComponent
                                            name={seqName}
                                            seq={sequence}
                                            viewer={viewerType}
                                            annotations={annotations}
                                            style={{ height: "100%", width: "100%" }}
                                        />
                                    ) : !SeqVizComponent && sequence ? (
                                        <div className="absolute inset-0 flex items-center justify-center text-xs text-muted-foreground">
                                            Loading sequence visualizer...
                                        </div>
                                    ) : (
                                        <div className="absolute inset-0 flex items-center justify-center text-xs text-muted-foreground">
                                            Please insert a sequence to visualize
                                        </div>
                                    )}
                                </div>
                            </div>

                            {/* Aligned Dot-Bracket Correlation Visualizer */}
                            <div className="bg-card text-card-foreground p-6 border border-border space-y-4">
                                <div className="flex items-center justify-between border-b border-border pb-3">
                                    <h2 className="font-bold text-card-foreground/80 text-xs tracking-wider uppercase flex items-center space-x-1.5">
                                        <Activity className="w-4 h-4 text-primary" />
                                        <span>Sequence & Dot-Bracket Correlation Map</span>
                                    </h2>
                                    <span className="text-[9px] text-muted-foreground font-extrabold uppercase px-2 py-0.5 bg-muted/60 border border-border">
                                        Interactive Aligner
                                    </span>
                                </div>

                                {/* Active Base Pairing status message */}
                                <div className="bg-muted/40 p-3 text-xs font-semibold flex items-center justify-between border border-border/40 min-h-[40px] px-4">
                                    {hoveredBaseIndex !== null ? (
                                        (() => {
                                            const baseChar = sequence[hoveredBaseIndex] || "";
                                            const dbChar = alignedStructure[hoveredBaseIndex] || "";
                                            const partnerIdx = basePairs[hoveredBaseIndex];
                                            const annotation = getAnnotationAt(hoveredBaseIndex);

                                            return (
                                                <div className="flex items-center justify-between w-full">
                                                    <div className="flex items-center space-x-2">
                                                        <span className="text-muted-foreground">Base:</span>
                                                        <span className="bg-primary/10 text-primary px-1.5 py-0.5 font-bold font-mono text-sm">{baseChar}</span>
                                                        <span className="text-muted-foreground">at position</span>
                                                        <span className="text-foreground font-bold">{hoveredBaseIndex + 1} bp</span>
                                                    </div>

                                                    <div className="flex items-center space-x-2">
                                                        {partnerIdx !== undefined ? (
                                                            <>
                                                                <span className="text-muted-foreground">Pairs with:</span>
                                                                <span className="bg-emerald-500/10 text-emerald-500 px-1.5 py-0.5 font-bold font-mono text-sm">
                                                                    {sequence[partnerIdx]}
                                                                </span>
                                                                <span className="text-muted-foreground">at</span>
                                                                <span className="text-emerald-500 font-bold font-mono text-sm">#{partnerIdx + 1} bp</span>
                                                                <span className="text-muted-foreground">({dbChar} ⟷ {alignedStructure[partnerIdx]})</span>
                                                            </>
                                                        ) : (
                                                            <span className="text-amber-500 font-bold bg-amber-500/10 px-2.5 py-0.5">
                                                                Unpaired Base (Loops/Bulge)
                                                            </span>
                                                        )}
                                                    </div>

                                                    {annotation && (
                                                        <div className="flex items-center space-x-1.5 border-l border-border/60 pl-3">
                                                            <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: annotation.color }} />
                                                            <span className="text-[10px] uppercase font-bold text-muted-foreground">{annotation.name}</span>
                                                        </div>
                                                    )}
                                                </div>
                                            );
                                        })()
                                    ) : (
                                        <span className="">
                                        </span>
                                    )}
                                </div>

                                {/* Aligning wrapped rows */}
                                <div className="space-y-6 overflow-x-auto py-2">
                                    {chunks.map((indices, chunkIdx) => (
                                        <div key={chunkIdx} className="space-y-1 min-w-[550px] border-b border-border/20 pb-4 last:border-0">
                                            {/* Column Index Markers */}
                                            <div className="flex text-[9px] font-bold text-muted-foreground font-mono">
                                                <div className="w-16 shrink-0 flex items-center justify-end pr-2 text-right uppercase">Position</div>
                                                {indices.map((idx) => (
                                                    <div key={`marker-${idx}`} className="w-7 text-center">
                                                        {(idx + 1) % 5 === 0 ? idx + 1 : ""}
                                                    </div>
                                                ))}
                                            </div>

                                            {/* Sequence Row */}
                                            <div className="flex items-center">
                                                <div className="w-16 shrink-0 text-[10px] font-bold text-muted-foreground uppercase tracking-wider pr-2 text-right">Sequence</div>
                                                <div className="flex">
                                                    {indices.map((idx) => {
                                                        const isHovered = idx === hoveredBaseIndex;
                                                        const isPartner = hoveredBaseIndex !== null && idx === basePairs[hoveredBaseIndex];
                                                        const ann = getAnnotationAt(idx);

                                                        let cellBg = "bg-muted/40 border-border/40 text-foreground";
                                                        if (isHovered) {
                                                            cellBg = basePairs[idx] !== undefined
                                                                ? "bg-primary/20 border-primary text-primary shadow-inner"
                                                                : "bg-amber-500/20 border-amber-500 text-amber-500 shadow-inner";
                                                        } else if (isPartner) {
                                                            cellBg = "bg-emerald-500/20 border-emerald-500 text-emerald-500";
                                                        }

                                                        return (
                                                            <div
                                                                key={`seq-${idx}`}
                                                                onMouseEnter={() => setHoveredBaseIndex(idx)}
                                                                onMouseLeave={() => setHoveredBaseIndex(null)}
                                                                className={`w-7 h-7 flex items-center justify-center font-mono text-xs font-bold border-r border-b cursor-pointer transition-all duration-150 ${cellBg}`}
                                                            >
                                                                {sequence[idx]}
                                                            </div>
                                                        );
                                                    })}
                                                </div>
                                            </div>

                                            {/* Structure Dot-Bracket Row */}
                                            <div className="flex items-center">
                                                <div className="w-16 shrink-0 text-[10px] font-bold text-muted-foreground uppercase tracking-wider pr-2 text-right">Structure</div>
                                                <div className="flex">
                                                    {indices.map((idx) => {
                                                        const isHovered = idx === hoveredBaseIndex;
                                                        const isPartner = hoveredBaseIndex !== null && idx === basePairs[hoveredBaseIndex];
                                                        const char = alignedStructure[idx];

                                                        let cellBg = "bg-muted/10 border-border/20 text-muted-foreground";
                                                        if (isHovered) {
                                                            cellBg = basePairs[idx] !== undefined
                                                                ? "bg-primary/10 border-primary text-primary"
                                                                : "bg-amber-500/10 border-amber-500 text-amber-500";
                                                        } else if (isPartner) {
                                                            cellBg = "bg-emerald-500/10 border-emerald-500 text-emerald-500";
                                                        }

                                                        return (
                                                            <div
                                                                key={`struct-${idx}`}
                                                                onMouseEnter={() => setHoveredBaseIndex(idx)}
                                                                onMouseLeave={() => setHoveredBaseIndex(null)}
                                                                className={`w-7 h-7 flex items-center justify-center font-mono text-xs font-extrabold border-r border-b cursor-pointer transition-all duration-150 ${cellBg}`}
                                                            >
                                                                {char}
                                                            </div>
                                                        );
                                                    })}
                                                </div>
                                            </div>

                                            {/* Feature Annotation Bands */}
                                            <div className="flex items-center">
                                                <div className="w-16 shrink-0 text-[9px] font-bold text-muted-foreground uppercase tracking-wider pr-2 text-right">Features</div>
                                                <div className="flex">
                                                    {indices.map((idx) => {
                                                        const ann = getAnnotationAt(idx);
                                                        return (
                                                            <div key={`annband-${idx}`} className="w-7 flex flex-col justify-start items-center">
                                                                {ann ? (
                                                                    <div
                                                                        className="h-1.5 w-full mt-0.5"
                                                                        style={{ backgroundColor: ann.color }}
                                                                        title={ann.name}
                                                                    />
                                                                ) : (
                                                                    <div className="h-1.5 w-full mt-0.5 bg-transparent" />
                                                                )}
                                                            </div>
                                                        );
                                                    })}
                                                </div>
                                            </div>

                                        </div>
                                    ))}
                                </div>
                            </div>

                            {/* Analytics & Computed Properties Card */}
                            <div className="bg-card text-card-foreground p-6 border border-border space-y-6">
                                <div className="flex items-center justify-between border-b border-border pb-3">
                                    <h2 className="font-bold text-card-foreground/80 text-xs tracking-wider uppercase flex items-center space-x-1.5">
                                        <Activity className="w-4 h-4 text-primary" />
                                        <span>Sequence Analytics & Computed Properties</span>
                                    </h2>
                                    <span className="text-[10px] text-muted-foreground font-semibold px-2 py-0.5 bg-muted/60 border border-border">
                                        Real-time Analysis
                                    </span>
                                </div>

                                <div className="grid grid-cols-1 md:grid-cols-3 gap-6">

                                    {/* Column 1: Core Biophysical Metrics */}
                                    <div className="space-y-4">
                                        <h3 className="text-[10px] font-bold text-muted-foreground uppercase tracking-wide border-b border-border/40 pb-1 flex items-center space-x-1">
                                            <Calculator className="w-3.5 h-3.5 text-primary/85" />
                                            <span>Biophysical Metrics</span>
                                        </h3>
                                        <div className="space-y-3">
                                            <div>
                                                <div className="flex justify-between text-[11px] font-semibold text-muted-foreground">
                                                    <span>GC Content</span>
                                                    <span className={gcContent >= 40 && gcContent <= 60 ? "text-emerald-500 font-bold" : "text-amber-500 font-bold"}>
                                                        {gcContent >= 40 && gcContent <= 60 ? "Optimal" : "Suboptimal"}
                                                    </span>
                                                </div>
                                                <div className="text-xl font-bold text-foreground mt-0.5">
                                                    {gcContent.toFixed(1)}%
                                                </div>
                                            </div>

                                            <div>
                                                <div className="flex justify-between text-[11px] font-semibold text-muted-foreground">
                                                    <span>Molecular Weight</span>
                                                    <span>{sequence.includes("U") ? "RNA" : "DNA"}</span>
                                                </div>
                                                <div className="text-xl font-bold text-foreground mt-0.5">
                                                    {(molWeight / 1000).toFixed(2)} kDa
                                                </div>
                                            </div>


                                        </div>
                                    </div>

                                    {/* Column 2: Nucleotide Frequencies */}
                                    <div className="space-y-4">
                                        <h3 className="text-[10px] font-bold text-muted-foreground uppercase tracking-wide border-b border-border/40 pb-1 flex items-center space-x-1">
                                            <Dna className="w-3.5 h-3.5 text-primary/85" />
                                            <span>Nucleotide Frequencies</span>
                                        </h3>
                                        <div className="space-y-2">
                                            {[
                                                { label: "Adenine (A)", count: freqs.A, pct: freqs.A_pct, color: "bg-red-500/80" },
                                                { label: sequence.includes("U") ? "Uracil (U)" : "Thymine (T)", count: freqs.T, pct: freqs.T_pct, color: "bg-amber-500/80" },
                                                { label: "Guanine (G)", count: freqs.G, pct: freqs.G_pct, color: "bg-emerald-500/80" },
                                                { label: "Cytosine (C)", count: freqs.C, pct: freqs.C_pct, color: "bg-blue-500/80" }
                                            ].map((item) => (
                                                <div key={item.label} className="text-xs space-y-1">
                                                    <div className="flex justify-between font-semibold">
                                                        <span className="text-muted-foreground">{item.label}</span>
                                                        <span className="text-foreground">{item.count} <span className="text-[10px] text-muted-foreground font-normal">({item.pct.toFixed(0)}%)</span></span>
                                                    </div>
                                                    <div className="w-full bg-muted/60 h-1.5 overflow-hidden">
                                                        <div className={`h-full ${item.color}`} style={{ width: `${item.pct}%` }} />
                                                    </div>
                                                </div>
                                            ))}
                                        </div>
                                    </div>

                                    {/* Column 3: Restriction Endonuclease Sites */}
                                    <div className="space-y-4">
                                        <h3 className="text-[10px] font-bold text-muted-foreground uppercase tracking-wide border-b border-border/40 pb-1 flex items-center space-x-1">
                                            <Scissors className="w-3.5 h-3.5 text-primary/85" />
                                            <span>Restriction Cut Sites</span>
                                        </h3>
                                        <div className="max-h-[140px] overflow-y-auto space-y-1.5 pr-1">
                                            {restrictionSites.length === 0 ? (
                                                <div className="text-[11px] text-muted-foreground text-center py-6 font-medium">
                                                    No restriction sites found (EcoRI, BamHI, etc.)
                                                </div>
                                            ) : (
                                                restrictionSites.map((site) => (
                                                    <div key={site.name} className="flex justify-between items-start bg-muted/40 p-1.5 border border-border/40 text-xs">
                                                        <div>
                                                            <div className="font-bold text-foreground">{site.name}</div>
                                                            <div className="text-[9px] font-semibold text-muted-foreground font-mono">{site.site}</div>
                                                        </div>
                                                        <div className="text-[10px] font-bold text-primary text-right max-w-[100px] break-words">
                                                            {site.positions.map(pos => `${pos}bp`).join(", ")}
                                                        </div>
                                                    </div>
                                                ))
                                            )}
                                        </div>
                                    </div>

                                </div>

                                {/* Bottom Wide Card: Translation */}
                                <div className="border-t border-border/60 pt-4 space-y-3">
                                    <h3 className="text-[10px] font-bold text-muted-foreground uppercase tracking-wide flex items-center space-x-1">
                                        <FileText className="w-3.5 h-3.5 text-primary/85" />
                                        <span>Translated Amino Acid Sequence (Translation Product)</span>
                                    </h3>
                                    <div className="bg-muted/40 p-3.5 border border-border font-mono text-[11px] leading-relaxed text-muted-foreground select-all max-h-[120px] overflow-y-auto break-all">
                                        {translation ? (
                                            <div className="space-y-1">
                                                <div>{translation}</div>
                                                <div className="text-[9px] text-muted-foreground font-semibold">
                                                    Total Residues: {translation.replace(/\*/g, "").length} aa {translation.endsWith("*") ? "(Stop Codon Terminated)" : ""}
                                                </div>
                                            </div>
                                        ) : (
                                            <span className="text-xs text-muted-foreground font-medium">
                                                No coding sequences found or sequence too short for translation.
                                            </span>
                                        )}
                                    </div>
                                </div>

                            </div>

                        </div>

                    </div>
                </div>

            </div>
        </div>
    );
}
