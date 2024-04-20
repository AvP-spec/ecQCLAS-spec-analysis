import os
path = os.path.dirname(os.path.abspath(__file__))
print("avp-las01 __init__.py writing package path in set_dir.txt")
print(path)
# print(os.listdir(path))
with open("set_dir.txt", "wt") as f:
    f.write(path)

