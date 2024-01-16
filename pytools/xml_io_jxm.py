import numpy as np
from io import StringIO
from xml.etree import cElementTree
from pandas import read_csv
import copy
import sys

def wrap(p, anchor, box):
    return p - box * np.round((p - anchor) / box)


def wrap_circular(p, box):
    s = np.sin(p / box * 2.0 * np.pi).mean(axis=0)
    c = np.cos(p / box * 2.0 * np.pi).mean(axis=0)
    anchor = np.arctan(s/c) / np.pi / 2.0 * box
    return p - box * np.round((p - anchor) / box)


class uxml(object):
    class read(object):
        '''
        usage: 

        subclass 'read':
            a = uxml.read('init.xml')
            natoms = a.natoms
            time_step = a.time_step
            pos = a.nodes['position'] # np.ndarray, float64, (natoms, 3)
            types = a.nodes['types']  # np.ndarray, string, (natoms,)
            image = a.nodes['image']  # np.ndarray, int64, (natoms, 3)
        '''


        def __init__(self, filename, needed=[]):
            tree = cElementTree.ElementTree(file=filename)
            root = tree.getroot()
            c = root[0]
            self.nodes = {}
            self.natoms = int(c.attrib['natoms'])
            self.time_step = int(c.attrib['time_step'])
            self.dimensions = int(c.attrib['dimensions'])

            # init nbonds, nangle, ndihs
            self.nbonds = 0
            self.nangles = 0
            self.ndihedrals = 0

            for e in c:
                if e.tag == 'box':
                    self.box = np.array([float(e.attrib['lx']), float(e.attrib['ly']), float(e.attrib['lz'])])
                    continue
                if ((len(needed) != 0) and (e.tag not in needed)):
                    continue
                if 'num' in e.attrib and e.attrib['num'] == '0':
                    continue

                if sys.version > '3':
                    self.nodes[e.tag] = read_csv(StringIO(e.text), delim_whitespace=True, header=None).squeeze("columns").values
                else:
                    self.nodes[e.tag] = read_csv(StringIO(unicode(e.text)), delim_whitespace=True, header=None).squeeze("columns").values

                if e.tag == 'bond':
                    self.nbonds = len(self.nodes['bond'])
                if e.tag == 'angle':
                    self.nangles = len(self.nodes['angle'])
                if e.tag == 'dihedral':
                    self.ndihedrals = len(self.nodes['dihedral'])
            
    class dump(object):
        '''

        subclass 'dump':
            User can load uxml.read class instance, make some revisions, and dump xml:
            ```
            o = uxml.dump(load=a)variab 'a' is a uxml.read class instance
            o.nodes['position'][:10, 0] = 0 # set particle 0-9, with position x = 0
            o.write('dump.xml')
            ```

            or this class instant can be built from nothing, and user input all info:

            ```
            o = uxml.dump()
            o.nodes['position'] = np.asarray([[0,0,0],[1,1,1],[2,2,1]], dtype=float)
            o.nodes['type'] = ['A', 'A', 'B']
            o.nodes['bond'] = [['A-A', 0, 1], ['A-B'], 1, 2]
            o.nodes['angle'] = [['polymer', 0, 1, 2]]
            ```

            then write to xml:

            ```
            o.write(outname='out.xml')
            ```

            some header info can be specified:

            ```
            o.write(outname='dump.xml', xml_type='galamost', allow=['position', 'type', 'bond'], time_step=0, box=[30, 30, 30]):
            ```
            params: outname  dump to xml named `outname`. Default, 'dump.xml'
            params: xml_type change the entry in xml for galamost/hoomd. Default, 'galamost'.
            params: allow  if allow is not specified (default as []), all the info in o.nodes dict will be dumped, or sub-items can be listd to dump. Default, [], all info.
            params: box  specified box, note that  if the dump instance come from nothing, `box` must be set, default None
            parmms: time_step  set time_step. Default 0

        '''

        def __init__(self, load=None):
            if load is not None:
                self.box = copy.deepcopy(load.box)
                self.time_step = copy.deepcopy(load.time_step)
                self.nodes = copy.deepcopy(load.nodes) 
            else:
                self.nodes = {}

        def write(self, outname='dump.xml', xml_type='galamost', allow=[], time_step=0, box=None):

            # check box 
            if not hasattr(self, 'box'):
                assert box != None, "No box info! set with write(box=[30, 30, 30])"
                self.box = box

            # set headers
            self.time_step = time_step
            self.natoms = len(self.nodes['position'])
            self.dimensions = 3 

            self.o = open(outname, 'w')
            self.formatter = {
                'position': '{:>18.8f}{:>18.8f}{:>18.8f}\n',
                'velocity': '{:>18.8f}{:>18.8f}{:>18.8f}\n',
                'bond': '{:<s} {:<d} {:<d}\n',
                'angle': '{:<s} {:<d} {:<d} {:<d}\n',
                'image': '{:<d} {:<d} {:<d}\n',
                'type': '{:<s}\n',
                'diameter': '{:<.4f}\n',
                'body': '{:<d}\n',
                'h_init': '{:<d}\n',
                'h_cris': '{:<d}\n',
                'mass': '{:<.4f}\n',
                'charge': '{:<.8f}\n',
                'dihedral': '{:<s} {:<d} {:<d} {:<d} {:<d}\n'
            }

            # 'write headers'
            self.o.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            self.o.write('<%s_xml version="1.6">\n' % (xml_type))
            self.o.write('<configuration time_step="%s" dimensions="%s" natoms="%s">\n' % (self.time_step, self.dimensions, self.natoms))
            self.o.write('<box lx="%f" ly="%f" lz="%f" xy="0" xz="0" yz="0"/>\n' % (self.box[0], self.box[1], self.box[2]))

            for node, value in self.nodes.items():
                if len(allow) != 0:
                    if node not in allow:
                        continue

                if len(value) != 0:
                    num = len(value)
                else:
                    continue
                    
                # write header
                self.o.write('<%s num="%d">\n' % (node, num))

                # write contents
                for p in value:
                    if node in ['position', 'image', 'bond', 'angle', 'dihedral', 'velocity']:
                        self.o.write(self.formatter[node].format(*p))
                    elif node in ['type', 'body', 'h_init', 'h_cris', 'mass', 'diameter', 'charge']:
                        self.o.write(self.formatter[node].format(p))
                    else:
                        print("oops! '%s' is not in the current formatter"%(node))

                # write tail
                self.o.write('</%s>\n' % (node))

            # 'write tails'
            self.o.write('</configuration>\n')
            self.o.write('</%s_xml>\n' % (xml_type))
            self.o.close()
            

    class slice(object):

        '''
        usasge: slice subsystem from xml object. Note that the sliced bond/angle/dihdral will be removed, and the index will be resorted.
        '''

        def __init__(self, xmlobject, indexs=[], types=[]):
            self.indexs = indexs
            self.u = xmlobject 
            self.node_particle = ['position', 'type', 'mass', 'image', 'charge', 'h_init', 'h_cris', 'diameter', 'body']
            self.node_topol = ['bond', 'angle', 'dihedral']

            # get the indexs
            self.id = []
            if len(indexs) != 0:
                self.id.extend(indexs)

            if len(types) != 0:
                index_types = np.argwhere(np.isin(self.u.nodes['type'], types)).reshape(-1)
                self.id.extend(index_types.tolist())

            self.id = np.unique(np.array(self.id, dtype=int))
            assert len(self.id) > 0,"oops! You do not specify any particles to be filtered!"

            self.u.natoms = len(self.id)

        def excute(self):
            for node in self.u.nodes:
                if node in self.node_particle:
                    self.u.nodes[node] = self.u.nodes[node][self.id]
                elif node in self.node_topol:

                    # get the elements contains self.id
                    index = np.unique(np.argwhere(np.isin(self.u.nodes[node][:, 1:], self.id))[:, 0])
                    inode = self.u.nodes[node][index]
                    row, col = inode.shape 

                    # remove bonds that contains atom not in self.id 
                    remove = np.isin(inode[:, 1:], self.id).sum(axis=-1)
                    inode = inode[remove == col-1]

                    # map the index in inode (bond, angle, ...), to new index in self.id
                    inode_index = inode[:, 1:].reshape(-1)
                    inode_index_new = np.argwhere(self.id == inode_index[:, None])[:, 1].reshape(-1, col-1)

                    inode[:, 1:] = inode_index_new
                    
                    # change the node in self.u
                    self.u.nodes[node] = inode
            return copy.deepcopy(self.u) 

if __name__ == '__main__':
    # read xml file into uxml.read intance
    xml = uxml.read('tests/a.xml')
    xml.nodes['position'][:, 0] = 0

    # load uxml.read instance to dump
    o = uxml.dump(load=xml)
    o.write("tests/out.xml")

    # build uxml.dump instance from nothing
    m = uxml.dump()
    m.box = [30, 30, 30]
    m.nodes['position'] = np.asarray([[0,0,0],[1,1,1],[2,2,1]], dtype=float)
    m.nodes['type'] = ['A', 'A', 'B']
    m.nodes['bond'] = np.asarray([['A-A', 0, 1], ['A-B', 1, 2]], dtype=object)
    m.nodes['angle'] = np.asarray([['polymer', 0, 1, 2]], dtype=object)
    m.write(outname='tests/build.xml', xml_type='galamost', allow=['position', 'type', 'bond'], time_step=1000, box=[40, 40, 40])

    # slice subsystem from xml instance 
    newxml = uxml.slice(xml, indexs = [], types=['A']).excute()
    uxml.dump(load=newxml).write('tests/slice.xml')
