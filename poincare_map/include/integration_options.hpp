#ifndef INTEGRATION_OPTIONS_HPP_6YVUVZRC
#define INTEGRATION_OPTIONS_HPP_6YVUVZRC


struct IntegrationOptions

{
    double abs_err = 1e-10;
    double rel_err = 1e-10;
    double dt = 1e-2;
    double coefficient = 1e3; //to be used with closing orbit detection strategies

    //Seemingly redundant. For usage with python bindings
    IntegrationOptions (double abs_error,
                        double rel_error,
                        double init_time_step,
                        double tol_coefficient);

    //Seemingly redundant. For usage with python bindings
    IntegrationOptions (double abs_error,
                        double rel_error,
                        double init_time_step);


    //Seemingly redundant. For usage with python bindings
    IntegrationOptions (double abs_error,
                        double rel_error);


    //Seemingly redundant. For usage with python bindings
    IntegrationOptions()=default;
};

#endif /* end of include guard: INTEGRATION_OPTIONS_HPP_6YVUVZRC */
