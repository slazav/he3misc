function check_dstr(dstr)
    assert(length(regexp(dstr, '^20\d\d\d\d\d\d$'))==1, '\nBad dstr: %s\n', dstr);
end
