def create_pwd_file(username, password, fname=None):
    if fname is None:
        fname = "{}_ieeglogin.bin".format(username[:3])
    with open(fname, "wb") as f:
        f.write(password.encode())
    print("-- -- IEEG password file saved -- --\n")
