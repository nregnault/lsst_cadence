from main_func import parser_args_set, main_set

argsp = parser_args_set()

main_set(argsp.folder, argsp.cadence)
