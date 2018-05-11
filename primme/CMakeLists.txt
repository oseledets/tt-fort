set(PRIMME_CSOURCES Complex.c
                    common_numerical.c
                    convergence_d.c
                    convergence_z.c
                    correction_d.c
                    correction_z.c
                    errors.c
                    factorize_d.c
                    factorize_z.c
                    init_d.c
                    init_z.c
                    inner_solve_d.c
                    inner_solve_z.c
                    locking_d.c
                    locking_z.c
                    main_iter_d.c
                    main_iter_z.c
                    numerical_d.c
                    numerical_z.c
                    ortho_d.c
                    ortho_z.c
                    primme_d.c
                    primme_f77.c
                    primme_interface.c
                    primme_z.c
                    restart_d.c
                    restart_z.c
                    solve_H_d.c
                    solve_H_z.c
                    update_projection_d.c
                    update_projection_z.c
                    update_W_d.c
                    update_W_z.c
                    wtime.c)
set(PRIMME_FSOURCES zdotcsub.f)
set(PRIMME_SOURCES ${PRIMME_CSOURCES} ${PRIMME_FSOURCES})

add_library(primme-obj OBJECT ${PRIMME_SOURCES})
add_library(primme-static STATIC $<TARGET_OBJECTS:primme-obj>)
add_library(primme-dynamic SHARED $<TARGET_OBJECTS:primme-obj>)
add_library(primme ALIAS primme-static)  # backward compatibility
