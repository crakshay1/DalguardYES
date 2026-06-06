import { useEffect, useRef, useState } from "react";
import type { Location } from "react-router-dom";
import { cn } from "../lib/utils";

type Props = {
  mainContentId: string;
  location: Location;
};

interface HeadingItem {
  id: string;
  text: string;
  level: number;
}

function extractHeadingsFromDOM(element: HTMLElement | null): HeadingItem[] {
  if (!element) return [];
  const headings: HeadingItem[] = [];
  let idCounter = 0;

  element.querySelectorAll("h1, h2, h3, h4").forEach((el) => {
    let id = el.id;
    if (!id) {
      id = `minimap-heading-${idCounter++}`;
      el.id = id;
    }
    const level = parseInt(el.tagName[1], 10) - 1;
    const text = el.textContent?.trim() || "";

    headings.push({
      id,
      text,
      level,
    });
  });
  return headings;
}

const LEVEL_STYLES = [
  "text-xs font-bold pl-2 text-primary",                          // h1 - Level 0
  "text-[11px] font-semibold pl-4 text-foreground/85",            // h2 - Level 1
  "text-[10px] font-medium pl-6 text-foreground/70",             // h3 - Level 2
  "text-[9px] pl-8 text-foreground/50",                           // h4 - Level 3
];

export function TextMinimap({ mainContentId, location }: Props) {
  const minimapRef = useRef<HTMLDivElement>(null);
  const [headings, setHeadings] = useState<HeadingItem[]>([]);
  const [activeId, setActiveId] = useState<string | null>(null);
  const [isMinimapHovered, setIsMinimapHovered] = useState(false);

  useEffect(() => {
    const timeout = setTimeout(() => {
      const main = document.getElementById(mainContentId);
      setHeadings(extractHeadingsFromDOM(main));
    }, 100);
    return () => clearTimeout(timeout);
  }, [mainContentId, location]);

  useEffect(() => {
    const main = document.getElementById(mainContentId);
    if (!main || !headings.length) return;

    const onScroll = () => {
      let found: string | null = null;
      for (let i = headings.length - 1; i >= 0; i--) {
        const h = document.getElementById(headings[i].id);
        if (h) {
          const rect = h.getBoundingClientRect();
          if (rect.top < window.innerHeight * 0.25) {
            found = headings[i].id;
            break;
          }
        }
      }
      setActiveId(found);

      if (found && minimapRef.current) {
        const el = minimapRef.current.querySelector(
          `[data-minimap-id='${found}']`,
        );
        if (el) {
          el.scrollIntoView({ block: "nearest" });
        }
      }
    };

    main.addEventListener("scroll", onScroll);
    window.addEventListener("scroll", onScroll);
    onScroll();

    return () => {
      main.removeEventListener("scroll", onScroll);
      window.removeEventListener("scroll", onScroll);
    };
  }, [mainContentId, headings]);

  const handleClick = (id: string) => {
    const el = document.getElementById(id);
    if (el) {
      const yOffset = -80;
      const y = el.getBoundingClientRect().top + window.pageYOffset + yOffset;
      window.scrollTo({ top: y, behavior: "smooth" });
    }
  };

  return (
    <div
      className="fixed right-6 top-28 z-30 font-mono bg-card/95 backdrop-blur-sm w-64 h-fit max-h-[75vh] overflow-y-auto border border-foreground/10 p-4 select-none cursor-pointer transition-all duration-200 rounded-none shadow-sm hover:shadow-md hover:border-foreground/20 hidden lg:block "
      ref={minimapRef}
      tabIndex={-1}
      aria-label="Minimap"
    >
      <div className="flex flex-col gap-1">
        {headings.map((h) => {
          const isActive = h.id === activeId;
          const levelStyle = LEVEL_STYLES[h.level] || LEVEL_STYLES[LEVEL_STYLES.length - 1];

          return (
            <div
              key={h.id}
              data-minimap-id={h.id}
              onClick={() => handleClick(h.id)}
              className={cn(
                "border-l-2 transition-all duration-150 py-1 pr-2 truncate rounded-none select-none",
                levelStyle,
                isActive
                  ? "bg-primary/10 border-primary font-bold text-primary dark:text-foreground"
                  : "border-transparent hover:bg-foreground/5 hover:border-foreground/15 hover:text-foreground dark:text-foreground"
              )}
              title={h.text}
            >
              {h.text}
            </div>
          );
        })}
      </div>
    </div>
  );
}

