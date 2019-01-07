#!/usr/bin/env python
# *-* coding: UTF-8 *-*

# Copyright 2012-2019 Ronald Römer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# export LD_LIBRARY_PATH=/home/zippy/VTK8/lib

import sys
sys.path.extend(['/home/zippy/VTK8/lib/python3.6/site-packages',
    '/home/zippy/vtkbool/build'])

import vtkboolPython
import vtk

import math
import os
import re

def extrude (pts, z, h):
    cell = vtk.vtkIdList()

    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetDataTypeToDouble()

    [ (vtk_pts.InsertNextPoint(pt[0], pt[1], z), cell.InsertNextId(i)) for i, pt in enumerate(pts) ]

    pd = vtk.vtkPolyData()
    pd.Allocate(1, 1)
    pd.SetPoints(vtk_pts)
    pd.InsertNextCell(vtk.VTK_POLYGON, cell)

    prod = vtk.vtkTrivialProducer()
    prod.SetOutput(pd)

    extr = vtk.vtkLinearExtrusionFilter()
    extr.SetInputConnection(prod.GetOutputPort())
    extr.SetVector(0, 0, h)

    pn = vtk.vtkPolyDataNormals()
    pn.SetInputConnection(extr.GetOutputPort())
    pn.AutoOrientNormalsOn()

    return pn

class Frieze:
    def __init__ (self, usr_cfg):
        self.cfg = { 'w': 39, 'a': 3.25/6, 'b': 1, 'c': .4333333333, 'e': .15/2, 'f': 3.25/3,
            'q': 1.0833333, 'r': .625, 'end': 'A', 'flip': False, 'ang': 45, 'end_seqs': [(), ()], 'shift': 0, 'sym': False }

        self.cfg.update(usr_cfg)

        if self.cfg['end'] == 'A':
            self.cfg.update({ 'o_a': self.cfg['a'], 'o_b': 0, 'o_c': 0 })

        elif self.cfg['end'] == 'B':
            self.cfg.update({ 'o_a': -self.cfg['b'], 'o_b': 0, 'o_c': self.cfg['b'] })

        elif self.cfg['end'] == 'C':
            self.cfg.update({ 'o_a': 0, 'o_b': self.cfg['a'], 'o_c': 0 })

        elif self.cfg['end'] == 'D':
            self.cfg.update({ 'o_a': -self.cfg['a'], 'o_b': 2*self.cfg['a'], 'o_c': 0 })

    def draw_bricks (self, seq, end_seq):
        seq = list(seq)

        pts = []

        if self.cfg['sym']:

            assert self.cfg['end'] == 'D'

            fixed = sum(seq[1:-1])+sum(end_seq)

            print(fixed)

            w = self.cfg['w']+2*self.cfg['a']

            i = 0

            while True:
                curr = (i*seq[-1]+fixed)*self.cfg['f']

                if (w-curr)/2 < seq[0]*self.cfg['f']:
                    break

                i += 1

            offs = (w-curr)/2

            seq_ = [0] + seq[1:-1] + [seq[-1]]*i + list(reversed(end_seq))

            for i in range(1, len(seq_)):
                seq_[i] += seq_[i-1]

            print(seq_)

            for s in seq_:
                mid = self.cfg['a']-offs-s*self.cfg['f']

                pts.extend([[mid+self.cfg['e'], -self.cfg['a']],
                    [mid+self.cfg['e'], -self.cfg['a']+2*self.cfg['e']],
                    [mid-self.cfg['e'], -self.cfg['a']+2*self.cfg['e']],
                    [mid-self.cfg['e'], -self.cfg['a']]])

        else:

            n, m = divmod(self.cfg['w']+self.cfg['o_b'], self.cfg['f'])

            if abs(m-self.cfg['f']) < 1e-5:
                n += 1
                m = 0

            print(n, m)

            n = int(n)
            n, m_ = divmod(n-sum(seq[:-1]), seq[-1])

            print('->', n, m_)

            seq.extend([seq[-1]]*(n-1))

            print(seq)

            if m_ == 0 and m < 1e-5:
                del seq[-1]

            if end_seq:
                seq[-len(end_seq):] = reversed(end_seq)

            # wenn 0, dann löschen
            seq = [ s for s in seq if s > 0 ]

            for i in range(1, len(seq)):
                seq[i] += seq[i-1]

            for s in seq:
                mid = self.cfg['a']-s*self.cfg['f']

                pts.extend([[mid+self.cfg['e'], -self.cfg['a']],
                    [mid+self.cfg['e'], -self.cfg['a']+2*self.cfg['e']],
                    [mid-self.cfg['e'], -self.cfg['a']+2*self.cfg['e']],
                    [mid-self.cfg['e'], -self.cfg['a']]])

        pts.extend([[-self.cfg['w']+self.cfg['o_a'], -self.cfg['a']],
            [-self.cfg['w']+self.cfg['o_a'], self.cfg['b']],
            [self.cfg['a'], self.cfg['b']],
            [self.cfg['a'], -self.cfg['a']]])

        return pts

    def draw_zz_bricks (self):
        j = int(self.cfg['shift']//self.cfg['c'])+1

        self.cfg['shift'] = self.cfg['shift']%self.cfg['c']

        n, m = divmod(self.cfg['w']+self.cfg['o_c']-self.cfg['shift'], self.cfg['c'])

        n = int(n)

        print(n, m)

        if self.cfg['sym']:
            self.cfg['shift'] = m/2

        pts = [ [-i*self.cfg['c']-self.cfg['shift'], ((i+j)%2)*self.cfg['c']] for i in range(n+1) ]

        if m > 1e-3:
            f = -1 if pts[-1][1] > 1e-5 else 1

            pts.append([ pts[-1][0]-m, pts[-1][1]+f*m ])

        pts.extend([ [-self.cfg['w']-self.cfg['o_c'], self.cfg['b']],
            [0, self.cfg['b']] ])

        if self.cfg['shift'] > 1e-5:
            y = self.cfg['c']-self.cfg['shift'] if j%2 == 0 else self.cfg['shift']

            pts.append([0, y])

        return pts

    def draw_spacer (self):
        pts = [[-self.cfg['w']+self.cfg['o_a'], 2*self.cfg['e']],
            [-self.cfg['w']+self.cfg['o_a'], self.cfg['b']],
            [self.cfg['a'], self.cfg['b']],
            [self.cfg['a'], 2*self.cfg['e']]]

        if self.cfg['end'] == 'B':
            pts[:1] = [[-self.cfg['w']-2*self.cfg['e'], 2*self.cfg['e']],
                [-self.cfg['w']-2*self.cfg['e'], -self.cfg['a']],
                [-self.cfg['w']-self.cfg['b'], -self.cfg['a']]]

        return pts

    def export (self, name):
        extr = extrude(self.draw_bricks(self.cfg['seqs'][0], self.cfg['end_seqs'][0]), self.cfg['q']/2, self.cfg['r'])
        extr1 = extrude(self.draw_spacer(), self.cfg['q']/2+self.cfg['r'], 2*self.cfg['e'])

        extr2 = extrude(self.draw_bricks(self.cfg['seqs'][1], self.cfg['end_seqs'][1]), -self.cfg['q']/2, -self.cfg['r'])
        extr3 = extrude(self.draw_spacer(), -self.cfg['q']/2-self.cfg['r'], -2*self.cfg['e'])

        extr4 = extrude(self.draw_zz_bricks(), -self.cfg['q']/2, self.cfg['q'])

        # extr + extr1
        bf = vtkboolPython.vtkPolyDataBooleanFilter()
        bf.SetInputConnection(extr.GetOutputPort())
        bf.SetInputConnection(1, extr1.GetOutputPort())
        bf.DecPolysOff()

        # extr2 + extr3
        bf1 = vtkboolPython.vtkPolyDataBooleanFilter()
        bf1.SetInputConnection(extr2.GetOutputPort())
        bf1.SetInputConnection(1, extr3.GetOutputPort())
        bf1.DecPolysOff()

        app = vtk.vtkAppendPolyData()
        app.AddInputConnection(bf.GetOutputPort())
        app.AddInputConnection(bf1.GetOutputPort())

        bf2 = vtkboolPython.vtkPolyDataBooleanFilter()
        bf2.SetInputConnection(app.GetOutputPort())
        bf2.SetInputConnection(1, extr4.GetOutputPort())
        bf2.DecPolysOff()

        ang = self.cfg['ang']*math.pi/180

        v = [0, 5]

        _v = [math.cos(ang)*v[0]-math.sin(ang)*v[1],
            math.sin(ang)*v[0]+math.cos(ang)*v[1]]

        plane = vtk.vtkPlaneSource()
        plane.SetOrigin(0, 0, 0)
        plane.SetPoint1(0, 0, 5)
        plane.SetPoint2(_v[0], _v[1], 0)
        plane.SetCenter(0, 0 , 0)

        bf3 = vtkboolPython.vtkPolyDataBooleanFilter()
        bf3.SetInputConnection(bf2.GetOutputPort())
        bf3.SetInputConnection(1, plane.GetOutputPort())
        bf3.SetOperModeToDifference()
        bf3.DecPolysOff()

        result = bf3

        if self.cfg['end'] == 'D':

            _v = [math.cos(-ang)*v[0]-math.sin(-ang)*v[1],
                math.sin(-ang)*v[0]+math.cos(-ang)*v[1]]

            plane1 = vtk.vtkPlaneSource()
            plane1.SetOrigin(0, 0, 0)
            plane1.SetPoint2(0, 0, 5)
            plane1.SetPoint1(_v[0], _v[1], 0)
            plane1.SetCenter(-self.cfg['w'], 0 , 0)

            bf4 = vtkboolPython.vtkPolyDataBooleanFilter()
            bf4.SetInputConnection(result.GetOutputPort())
            bf4.SetInputConnection(1, plane1.GetOutputPort())
            bf4.SetOperModeToDifference()
            bf4.DecPolysOff()

            result = bf4

        if 'clip' in self.cfg:
            plane2 = vtk.vtkPlaneSource()
            plane2.SetOrigin(0, 0, 0)
            plane2.SetPoint1(0, 0, 5)
            plane2.SetPoint2(v[0], v[1], 0)
            plane2.SetCenter(-self.cfg['clip'], 0 , 0)

            bf5 = vtkboolPython.vtkPolyDataBooleanFilter()
            bf5.SetInputConnection(result.GetOutputPort())
            bf5.SetInputConnection(1, plane2.GetOutputPort())
            bf5.SetOperModeToDifference()
            bf5.DecPolysOff()

            result = bf5

        if 'pins' in self.cfg:
            app1 = vtk.vtkAppendPolyData()

            t = self.cfg.get('clip', 0)
            fake_w = self.cfg.get('fake_w', self.cfg['w'])

            u = (fake_w-t)/self.cfg['pins']/2

            for i in range(self.cfg['pins']):
                mid = t+u*(1+i*2)

                pin = extrude([ [-mid-2.5, self.cfg['b']], [-mid-2.5, self.cfg['b']+1.5],
                    [-mid+2.5, self.cfg['b']+1.5], [-mid+2.5, self.cfg['b']] ], -.75, 1.5)

                app1.AddInputConnection(pin.GetOutputPort())

            bf6 = vtkboolPython.vtkPolyDataBooleanFilter()
            bf6.SetInputConnection(result.GetOutputPort())
            bf6.SetInputConnection(1, app1.GetOutputPort())
            bf6.DecPolysOff()

            result = bf6

        tra = vtk.vtkTransform()
        tra.Scale(-1, 1, 1)

        tf = vtk.vtkTransformPolyDataFilter()
        tf.SetInputConnection(result.GetOutputPort())
        tf.SetTransform(tra)

        if self.cfg['flip']:
            tf.Update()
            pd = tf.GetOutput()

            for i in range(pd.GetNumberOfCells()):
                pd.ReverseCell(i)

            result = tf

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(result.GetOutputPort())

        writer = vtk.vtkPolyDataWriter()
        writer.SetInputConnection(clean.GetOutputPort())
        writer.WriteToOutputStringOn()
        writer.Update()

        dat = writer.GetOutputString()

        pd = clean.GetOutput()
        num = pd.GetNumberOfPoints()

        pts = [ [ '%.8f' % x for x in pd.GetPoint(i) ] for i in range(num) ]

        m = re.search('POINTS.*?\n(.*?)\nP', dat, re.S)

        with open(name, 'w') as f:
            f.write(dat[:m.start(1)] + ' '.join(sum(pts, [])) + dat[m.end(1):])

        if os.path.exists('stl'):
            tri = vtk.vtkTriangleFilter()
            tri.SetInputConnection(clean.GetOutputPort())

            writer1 = vtk.vtkSTLWriter()
            writer1.SetFileName('stl/{0}.stl'.format(name[:-4]))
            writer1.SetInputConnection(tri.GetOutputPort())
            writer1.Update()

if __name__ == '__main__':
    cfgs = [
        { 'seqs': [(2, 2), (2, 1, 2)], 'w': 39., 'pins': 2, 'fake_w': 40*3.25/3 },
        { 'seqs': [(1, 1), (1, 1)], 'w': 3.25, 'end': 'B', 'flip': True },
        { 'seqs': [(2, 2), (2, 1, 2)], 'w': 40*3.25/3, 'pins': 2 },
        { 'seqs': [(1, 1), (1, 1)], 'w': 11.375, 'end': 'C', 'flip': True  },
        { 'seqs': [(2, 2), (2, 1, 2)], 'w': 91*3.25/3, 'end': 'D', 'end_seqs': [(), (1,)], 'pins': 5 },
        { 'seqs': [(1, 1), (1, 1)], 'w': 85*3.25/3, 'end': 'C', 'end_seqs': [(0,), (0,)], 'shift': 3*.4333333333/2, 'pins': 5 },
        { 'seqs': [(2, 1, 2), (2, 2)], 'w': 40*3.25/3, 'end': 'C', 'flip': True, 'end_seqs': [(), (1, 2)], 'pins': 2 },
        { 'seqs': [(1, 1), (1, 1)], 'w': 3.25, 'end': 'B' },
        { 'seqs': [(2, 2), (2, 1, 2)], 'w': 40*3.25/3, 'flip': True, 'clip': 8*3.25/3, 'pins': 2 },
        { 'seqs': [(1, 1), (1, 1)], 'w': 7*3.25/3, 'end': 'B', 'flip': True, 'ang': 22.5, 'pins': 1 },
        { 'seqs': [(1, 1), (1, 1)], 'w': 47*3.25/3, 'ang': 22.5, 'pins': 3 },
        { 'seqs': [(2, 2), (2, 1, 2)], 'w': 18.3848, 'end': 'D', 'ang': 22.5, 'sym': True, 'end_seqs': [(1,), ()] },
        { 'seqs': [(1,), (1,)], 'w': 3.25/3, 'end': 'B', 'flip': True },
        { 'seqs': [(1, 1), (1, 1)], 'w': 11.5*3.25/3, 'end': 'C', 'pins': 1 },
        { 'seqs': [(1, 1), (1, 1)], 'w': 11.5*3.25/3, 'end': 'C', 'flip': True, 'pins': 1 }
    ]

    for i, cfg in enumerate(cfgs):
        print('~~ {0} ~~'.format(i))

        Frieze(cfg).export('test{0}.vtk'.format(i))
