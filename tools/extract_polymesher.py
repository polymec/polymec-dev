def find_dependencies(source_file, deps = set()):
    import os.path
    f = open(source_file)
    lines = f.readlines()
    f.close()
    print 'searching %s'%source_file
    for line in lines:
        if '#include "' in line:
            header = line[10:-2]
            if os.path.exists(header) and header not in deps:
                deps.add(header)
                more_deps = find_dependencies(header, deps)
                for dep in more_deps:
                    deps.add(dep)
            source = header.replace('.h', '.c')
            if os.path.exists(source) and source not in deps:
                deps.add(source)
                more_deps = find_dependencies(source, deps)
                for dep in more_deps:
                    deps.add(dep)
    return deps

print find_dependencies('polymesher/polymesher.c')
