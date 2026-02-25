# import logging

# logging.basicConfig(
#     level=logging.INFO,
#     format="[%(levelname)s] %(message)s"
# )
# logger = logging.getLogger(__name__)



# def get_logger(name="interface_constructor", level=logging.INFO):

#     logger = logging.getLogger(name)
    
#     if not logger.handlers:
#         handler = logging.StreamHandler()
#         formatter = logging.Formatter(
#             "[%(levelname)s] %(message)s"
#         )
#         handler.setFormatter(formatter)
#         logger.addHandler(handler)
#     logger.setLevel(level)
    
#     return logger
