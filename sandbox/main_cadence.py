from main_func import main_cadence, parser_args_cadence

argsp = parser_args_cadence()

main_cadence(argsp.data, argsp.folder, argsp.cadence)

