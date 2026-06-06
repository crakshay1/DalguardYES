export default function Footer() {

  const scrollToTop = () => {
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  return (
    <footer className="bg-background border-t border-foreground/10 py-10 px-8">
      <div className="max-w-5xl mx-auto flex flex-col md:flex-row justify-between items-center gap-6">


        <div className="text-center md:text-left">
          <h3 className="text-sm font-medium text-foreground mb-1">Dalguard</h3>
          <p className="text-xs text-foreground/50">
            Developed by students at Université d'Évry Paris-Saclay
          </p>
        </div>

      </div>
    </footer>
  );
}