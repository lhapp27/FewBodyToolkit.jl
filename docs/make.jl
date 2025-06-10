using FewBodyToolkit, Documenter, Literate, DocumenterCitations, Plots, Antique, Interpolations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

# Convert example scripts to markdown
example_files = [
  "example1D.jl",
  "example2D.jl",
  "example3D.jl",
  "3B1D_23body.jl",
  "1D_2+1.jl",
  "ISGL_HD+.jl",
]

for examplename in example_files
  exname = examplename[1:end-3]  # Remove the ".jl" extension
  Literate.markdown(joinpath(@__DIR__, "..", "examples", examplename), joinpath(@__DIR__, "src"), name=exname, flavor=Literate.DocumenterFlavor())
end



makedocs(
  sitename   = "FewBodyToolkit.jl",
  modules    = [FewBodyToolkit],
  format     = Documenter.HTML(collapselevel=1), #mathengine=MathJax()
  pagesonly  = true,
  pages      = [
    "Home" => "index.md",
    "Basis Functions" => "BasisFunctions.md",
     "Examples" => [
            "Overview" => "examples.md",
            "2B: Example 1D" => "example1D.md",
            "2B: Example 2D" => "example2D.md",
            "2B: Example 3D" => "example3D.md",
            "3B1D: Consistency with 2-body" => "3B1D_23body.md",
            "3B1D: 2+1 system" => "1D_2+1.md",
            "3B3D: HD+ system" => "ISGL_HD+.md",
        ],
    "Advanced Options" => "AdvancedOptions.md",
    "API"  => "api.md",
    "References" => "references.md",
  ],
  plugins=[bib],
)

deploydocs(
    repo = "github.com/lhapp27/FewBodyToolkit.jl.git",
    devbranch = "main",
    push_preview = true,
    target = "build",
)