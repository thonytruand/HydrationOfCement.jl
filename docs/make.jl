using Documenter, HydrationOfCement
 
makedocs(modules=[HydrationOfCement],
        doctest=true)
 
deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/thonytruand/HydrationOfCement.git",
    julia  = "1.6.2",
    osname = "windows")