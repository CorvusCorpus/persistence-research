# Build the write-up (needs a TeX installation with latexmk, e.g. texlive-latex-recommended + latexmk).
chain_sieve.pdf: chain_sieve.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error chain_sieve.tex

clean:
	latexmk -C chain_sieve.tex

.PHONY: clean
