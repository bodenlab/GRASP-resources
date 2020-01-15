def readOBOFile(obofile):
    """
    Read/load OBO file that contains ontology defs for
    Uber anatomy ontology (Uberon), Cell Ontology (CL) and Experimental Factor Ontology (EFO)
    http://cellontology.org
    see also http://www.obofoundry.org/

    Example of one "term" and one "typedef" entry (note CL refers to cell ontology and UBERON to the anatomy ontology:

    [Term]
    id: CL:0000513
    name: cardiac muscle myoblast
    namespace: cell
    alt_id: CL:0000714
    def: "A precursor cell destined to differentiate into cardiac muscle cell." [GOC:tfm, MESH:A11.635.470]
    synonym: "cardiac muscle progenitor cell" EXACT []
    synonym: "cardiomyocyte progenitor cell" EXACT []
    xref: FMA:84797
    is_a: CL:0002494 ! cardiocyte
    is_a: CL:0010021 ! cardiac myoblast
    intersection_of: CL:0000056 ! myoblast
    intersection_of: develops_into CL:0000746 ! cardiac muscle cell
    intersection_of: part_of UBERON:0001133 ! cardiac muscle tissue
    relationship: develops_into CL:0000746 ! cardiac muscle cell
    relationship: part_of UBERON:0001133 ! cardiac muscle tissue

    [Typedef]
    id: part_of
    name: part of
    def: "a core relation that holds between a part and its whole" []
    xref: BFO:0000050
    is_transitive: true
    is_a: overlaps ! overlaps
    inverse_of: has_part ! has part
    """
    src = open(obofile, 'r')
    terms = {}
    in_term_def = False
    in_type_def = False
    for line in src:
        if in_term_def:
            word = line.split()
            if line.startswith('id: '):
                term_id = word[1].strip()
                term_is = set()
            elif line.startswith('name: '):
                term_name = line[6:].strip()
            elif line.startswith('def: '):
                # Note this is a multi-line field, delimited by "'s
                pass
            elif line.startswith('is_a: '):
                term_is.add((word[1].strip(), 'is_a'))
            elif line.startswith('relationship: '):
                term_is.add((word[2], word[1]))
            elif line.startswith('intersection_of: '):
                pass
            elif line.startswith('is_obsolete: '):
                in_term_def = False # ignore this entry
        if line.startswith('[Term]'):
            if in_term_def: # already defining one, stash it before moving on to the next...
                terms[term_id] = (term_name, term_is)
            elif in_type_def:
                in_type_def = False
            in_term_def = True
        if line.startswith('[Typedef]'):
            if in_term_def: # already defining one, stash it before moving on to the next...
                terms[term_id] = (term_name, term_is)
                in_term_def= False
            in_type_def = True
    if in_term_def: #  defining one, stash it
        terms[term_id] = (term_name, term_is)
    return terms

def getChildren(terms, parent):
    all = []
    for t in terms:
        (name, isa) = terms[t]
        for (id, rel) in isa:
            if id == parent:
                all.append((t, rel))
    return all

def getParents(terms, child):
    all = []
    (name, isa) = terms[child]
    for (id, rel) in isa:
        all.append((id, rel))
    return all

def getID(terms, query_name):
    for t in terms:
        (name, isa) = terms[t]
        if name == query_name:
            return t

def getName(terms, id):
    return terms[id][0]

def listParents(terms, query_name):
    child_ID = getID(terms, query_name)
    all = []
    for (parent, rel) in getParents(terms, child_ID):
        all.append(getName(terms, parent))
    return all

if __name__ == '__main__':
    terms = readOBOFile('/Users/mikael/simhome/share/cl.obo')
    print(len(terms))