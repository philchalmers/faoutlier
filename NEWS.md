Changes in version 0.7

- added a `progress` argument to `gCD()`, `LD()`, and `GOF()` to print the progress of the iterations

Changes in version 0.5

  - several orgainizational overahauls

  - split LD() into LD() and GOF() for better clarity

  - avoid masking lavaan::sem with sem::sem

Changes in version 0.4

  - better use of the parallel package for faster computing in most functions (most notably
  in the forward.search implementation)

  - numerical tests added for better package stability

  - more flexible passing of lavaan-type models and other functions used within the package