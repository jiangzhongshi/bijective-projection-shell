#ifndef PRISM_EXTRACTION_HPP
#define PRISM_EXTRACTION_HPP

class PrismCage;
namespace prism {
  [[deprecated]] bool mid_surface_extraction(PrismCage&);
  bool shell_extraction(PrismCage& pc, bool base);
}

#endif