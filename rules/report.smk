# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

rule report:
    message: "Compile report."
    params:
        fn="report"
    input:
        tex="report/report.tex",
        bib="report/references.bib"
    output: "report/report.pdf"
    shell:
        """
        pdflatex -output-directory report {input.tex}
        cd report; bibtex {params.fn}; cd ..
        pdflatex -output-directory report {input.tex}
        pdflatex -output-directory report {input.tex}
        """