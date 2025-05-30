using Documenter, FewBodyToolkit, Literate, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

# Convert example scripts to markdown
#Literate.markdown("../examples/example1D.jl", joinpath(@__DIR__, "src"), name="example1D",flavor=Literate.DocumenterFlavor())
#Literate.markdown("../examples/example2D.jl", joinpath(@__DIR__, "src"), name="example2D",flavor=Literate.DocumenterFlavor())
#Literate.markdown("../examples/example3D.jl", joinpath(@__DIR__, "src"), name="example3D",flavor=Literate.DocumenterFlavor())

makedocs(
  sitename   = "FewBodyToolkit.jl",
  modules    = [FewBodyToolkit],
  format     = Documenter.HTML(), #mathengine=MathJax()
  pagesonly  = true,
  pages      = [
    "Home" => "index.md",
    "Basis Functions" => "BasisFunctions.md",
    "Example 1D" => "example1D.md",
    #"Example 2D" => "example2D.md",
    #"Example 3D" => "example3D.md",
    "API"  => "api.md",
    "References" => "references.md",
  ],
  plugins=[bib],
)
