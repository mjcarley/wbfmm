(TeX-add-style-hook
 "notes"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "times"
    "amsmath")
   (TeX-add-symbols
    "D"
    "E"
    "J")
   (LaTeX-add-bibliographies
    "abbrev"
    "maths"))
 :latex)

