#!/usr/bin/env python


def prosite_to_re(pattern):
    _prosite_trans = str.maketrans("abcdefghijklmnopqrstuvwxyzX}()<>",
                                  "ABCDEFGHIJKLMNOPQRSTUVW.YZ.]{}^$", "-.")

    """convert a valid Prosite pattern into an re str"""
    
    flg = (pattern[:2] == "[<")
    s = str.replace(pattern, "{", "[^")
    s = s.translate(_prosite_trans)
    # special case "[<" and ">]", if they exist
    if flg:
        i = str.index(s, "]")
        s = "(?:^|[" + s[2:i] + "])" + s[i+1:]
    if s[-2:] == "$]":
        i = str.rindex(s, "[")
        s = s[:i] + "(?:" + s[i:-2] + "]|$)"
    elif s[-3:] == "$]$":
        i = str.rindex(s, "[")
        s = s[:i] + "(?:" + s[i:-3] + "]|$)$"
    return s



#x = prosite_to_re("H-x(2,4)-[SC]-x(2)-{A}-x-[LIVMF](2)-[ST]-H-G")

#print(x)
