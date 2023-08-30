## Masking of repetitive elements

The polished reference genome is masked for repetitive elements in order to prepare the computational annotation.

### Script order:

- code/repeat_masking/run_repeatmodeler.sh
run Repeatmodeler to identify repetitive elements in the polished reference genome

- code/repeat_masking/run_repeatmasker.sh
run Repeatmasker on the Repeatmodeler output to classify identified elements and mask the polished reference genome

- code/repeat_masking/analyse_repetitive_elements.Rmd
analyse the repetitive element composition
