#include "integration_options.hpp"


IntegrationOptions::IntegrationOptions (double abs_error,
                                        double rel_error,
                                        double init_time_step,
                                        double tol_coefficient)
        : abs_err{abs_error},
          rel_err{rel_error},
          dt{init_time_step},
          coefficient{tol_coefficient} { }

IntegrationOptions::IntegrationOptions (double abs_error,
                                        double rel_error,
                                        double init_time_step)
        : abs_err{abs_error},
          rel_err{rel_error},
          dt{init_time_step} { }

IntegrationOptions::IntegrationOptions (double abs_error,
                                        double rel_error)
        : abs_err{abs_error},
          rel_err{rel_error}{ }

