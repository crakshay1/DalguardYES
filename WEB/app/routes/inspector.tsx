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

export default function Inspector() {
    const navigate = useNavigate();
    const [seqName, setSeqName] = useState("RiboGuard_Design_01");
    const [SeqVizComponent, setSeqVizComponent] = useState<any>(null);

    useEffect(() => {
        // Dynamically load SeqViz client-side to prevent SSR failures
        import("seqviz")
            .then((mod) => {
                setSeqVizComponent(() => mod.SeqViz);
            })
            .catch((err) => {
                console.error("Failed to load seqviz library:", err);
            });
    }, []);
    const [sequence, setSequence] = useState(
        "ATGCTAGCGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    );
    const [viewerType, setViewerType] = useState<"linear" | "circular" | "both">("linear");

    // Custom annotations state
    const [annotations, setAnnotations] = useState<Annotation[]>([
        { name: "Start Codon", start: 0, end: 3, direction: 1, color: "#3b82f6" },
        { name: "RBS Site", start: 12, end: 20, direction: 1, color: "#ef4444" },
        { name: "Spacer", start: 20, end: 35, direction: 1, color: "#eab308" },
    ]);

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
                                className={`w-full flex items-center space-x-3 px-4 py-3 rounded-lg text-sm font-semibold transition-all duration-150 ${item.active
                                        ? "bg-sidebar-primary text-sidebar-primary-foreground shadow-sm"
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
                            Structure Inspector
                        </h1>
                        <p className="text-xs text-muted-foreground mt-1 font-medium">
                            Analyze custom DNA/RNA/Protein sequences using LatchBio SeqViz visualization.
                        </p>
                    </div>

                    <div className="flex items-center space-x-4">
                        <ThemeToggle />
                    </div>
                </header>

                {/* Workspace Grid */}
                <div className="flex-1 overflow-y-auto p-8 space-y-6 bg-muted/20">
                    <div className="grid grid-cols-12 gap-6 max-w-[1400px] mx-auto items-stretch h-full">

                        {/* Left Controls Card (5 columns) */}
                        <div className="col-span-5 bg-card text-card-foreground p-6 rounded-xl border border-border shadow-sm flex flex-col justify-between space-y-4">
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
                                            className="w-full bg-muted/40 border border-border rounded-lg px-3 py-2 text-sm font-semibold text-foreground focus:outline-none focus:ring-1 focus:ring-primary"
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
                                            className="w-full bg-muted/40 border border-border rounded-lg px-3 py-2 text-xs font-mono uppercase text-foreground focus:outline-none focus:ring-1 focus:ring-primary leading-normal resize-none"
                                        />
                                        <span className="text-[10px] text-muted-foreground font-semibold block text-right">
                                            Length: {sequence.length} bp
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
                                                    className={`py-1.5 px-3 rounded-lg text-xs font-bold capitalize transition-all border ${viewerType === opt
                                                            ? "bg-primary text-primary-foreground border-primary"
                                                            : "bg-card border-border hover:bg-muted/30"
                                                        }`}
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
                                <div className="max-h-[150px] overflow-y-auto space-y-1.5 mb-3 border border-border/40 p-2 rounded-lg bg-muted/10">
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

                        {/* Right SeqViz Viewer Card (7 columns) */}
                        <div className="col-span-7 bg-card text-card-foreground p-6 rounded-xl border border-border shadow-sm flex flex-col justify-between items-stretch min-h-[550px]">
                            <div className="flex items-center justify-between border-b border-border pb-3 mb-4 shrink-0">
                                <h2 className="font-bold text-card-foreground/80 text-xs tracking-wider uppercase">
                                    Sequence Visualization
                                </h2>
                                <span className="text-[10px] text-muted-foreground font-semibold px-2 py-0.5 bg-muted/60 border border-border rounded">
                                    SeqViz Renderer
                                </span>
                            </div>

                            <div className="flex-1 bg-white dark:bg-[#1a2024]/40 rounded-xl overflow-hidden border border-border p-4 relative min-h-[480px]">
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

                    </div>
                </div>

            </div>
        </div>
    );
}
