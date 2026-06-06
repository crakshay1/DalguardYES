export interface FastaResult {
  name: string;
  sequence: string;
}

export interface CSVParseResult {
  headers: string[];
  rows: Record<string, string>[];
}

export interface CandidateData {
  rbs: string;
  spacer: string;
  orthScore: string;
  wtLeakage: string;
  rbsAccess: string;
  structure: string;
}

export interface ScatterDataPoint {
  id: number;
  wtLeakage: number;
  binding: number;
  access: number;
}

export interface DashboardInputData {
  orthogonalAntiSD?: string;
  wtAntiSD?: string;
  cdsStart?: string;
  targetExpression?: string;
  antiSD?: string;
  beforeRBS?: string;
  afterRBS?: string;
}

export interface DashboardJSONResult {
  inputs?: DashboardInputData;
  candidates?: CandidateData[];
  scatterPoints?: ScatterDataPoint[];
}

export type CSVIngestionType =
  | { type: "candidates"; data: CandidateData[] }
  | { type: "scatter"; data: ScatterDataPoint[] }
  | { type: "unknown"; headers: string[] };

/**
 * Standard CSV Parser helper
 */
export function parseCSV(text: string): CSVParseResult | null {
  const lines = text.split("\n").map(l => l.trim()).filter(Boolean);
  if (lines.length < 2) return null;

  const headers = lines[0].split(",").map(h => h.trim().toLowerCase());
  const rows: Record<string, string>[] = [];

  for (let i = 1; i < lines.length; i++) {
    const values = lines[i].split(",").map(v => v.trim());
    const row: Record<string, string> = {};
    headers.forEach((header, index) => {
      row[header] = values[index] || "";
    });
    rows.push(row);
  }

  return { headers, rows };
}

/**
 * Parses FASTA formatted text
 */
export function parseFASTA(text: string, defaultName: string): FastaResult {
  if (text.startsWith(">")) {
    const lines = text.split("\n");
    const header = lines[0].substring(1).trim();
    const sequence = lines.slice(1).join("").replace(/[^a-zA-Z]/g, "").toUpperCase();
    return { name: header || defaultName, sequence };
  } else {
    const sequence = text.replace(/[^a-zA-Z]/g, "").toUpperCase();
    return { name: defaultName, sequence };
  }
}

/**
 * Parses and classifies CSV contents based on headers
 */
export function ingestCSV(text: string): CSVIngestionType {
  const parsed = parseCSV(text);
  if (!parsed) {
    return { type: "unknown", headers: [] };
  }

  const { headers, rows } = parsed;

  if (headers.includes("rbs") && headers.includes("spacer")) {
    const data: CandidateData[] = rows.map((row) => ({
      rbs: row.rbs || "",
      spacer: row.spacer || "",
      orthScore: row.orthscore || row["orth score"] || "0.0",
      wtLeakage: row.wtleakage || row["wt leakage"] || "0.0",
      rbsAccess: row.rbsaccess || row["rbs access"] || "0.0",
      structure: row.structure || ".".repeat((row.rbs || "").length + (row.spacer || "").length + 10),
    }));
    return { type: "candidates", data };
  }

  if (headers.includes("wtleakage") || headers.includes("wt leakage") || headers.includes("binding") || headers.includes("orthogonal binding")) {
    const data: ScatterDataPoint[] = rows.map((row, idx) => {
      const wtLeakage = parseFloat(row.wtleakage) || parseFloat(row["wt leakage"]) || 0.0001;
      const binding = parseFloat(row.binding) || parseFloat(row["orthogonal binding"]) || 0.0;
      const access = parseFloat(row.access) || parseFloat(row.accessibility) || 0.5;
      return {
        id: idx,
        wtLeakage,
        binding,
        access,
      };
    });
    return { type: "scatter", data };
  }

  return { type: "unknown", headers };
}

/**
 * Standard JSON Ingest Parser
 */
export function ingestJSON(text: string): DashboardJSONResult {
  try {
    return JSON.parse(text) as DashboardJSONResult;
  } catch (err) {
    throw new Error("Failed to parse JSON. Ensure file is correctly formatted.");
  }
}
