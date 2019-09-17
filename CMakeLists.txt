cmake_minimum_required(VERSION 3.15)
project(ode_integrators VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)


# add_subdirectory(stand_alone_odeint_examples)



############## poincare_map ###############


include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/poincare_map)


add_subdirectory(poincare_map)


# Install the export set for use with the install-tree
install(EXPORT poincare_mapTargets
        FILE
        poincare_mapTargets.cmake
        NAMESPACE
        poincare_map::
        DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/cmake/poincare_map)


#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/poincare_mapConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)


configure_package_config_file (cmake/poincare_mapConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/poincare_mapConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR} )


#Install the config and configversion files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/poincare_mapConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/poincare_mapConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )



##############################################
## Exporting from the build tree


export(EXPORT poincare_mapTargets FILE ${CMAKE_CURRENT_BINARY_DIR}/poincare_mapTargets.cmake NAMESPACE poincare_map::)

#Register package in user's package registry
export(PACKAGE poincare_map)