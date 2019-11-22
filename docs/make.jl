using Documenter, ClinicalTrialUtilities

makedocs(
    modules = [ClinicalTrialUtilities],
    sitename = "ClinicalTrialUtilities",
    authors = "Vladimir Arnautov",
    linkcheck = false,
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Sample size" => "samplesize.md",
        "Power" => "power.md",
        "Confidence intervals" => "ci.md",
        "Descriptive statistics" => "ds.md",
        "NCA" => ["NCA" => "nca.md",
        "Pharmakokinetics" => "pk.md",
        "Pharmacodynamics" => "pd.md"],
        "Randomization" => "random.md",
        "Simulations" => "sim.md",
        "Utilities" => "utils.md",
        "Examples" => "examples.md",
        "Export" => "export.md",
        "Validation" => "validation.md",
        "References" => "ref.jl",
    ],
)

deploydocs(repo = "github.com/PharmCat/ClinicalTrialUtilities.jl.git")
