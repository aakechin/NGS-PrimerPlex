from kivy.app import App
from kivy.uix.button import Button
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.screenmanager import ScreenManager,Screen
from kivy.uix.gridlayout import GridLayout
from kivy.uix.floatlayout import FloatLayout
from kivy.core.window import Window
from kivy.config import ConfigParser
from kivy.uix.textinput import TextInput
from kivy.uix.label import Label
from kivy.uix.popup import Popup
from kivy.uix.filechooser import FileChooserListView
from kivy.uix.progressbar import ProgressBar
from kivy.metrics import dp
from kivy.graphics import Color,Rectangle
from kivy.properties import ObjectProperty
from kivy.properties import NumericProperty
from kivy.lang import Builder
from kivy.config import Config
from kivy.clock import Clock
from datetime import datetime
import os
import re
import time
import docker
import threading
import json
from docker.types import Mount
from functools import partial

Config.set('input', 'mouse', 'mouse,multitouch_on_demand')
global_image_names=['aakechin/ngs-primerplex:native',
                    'aakechin/ngs-primerplex:latest',
                    'ngs_primerplex_ref:latest']
available_images=[False,False,False]

class LoadDialog(FloatLayout):
    load = ObjectProperty(None)
    cancel = ObjectProperty(None)

class SaveDialog(FloatLayout):
    save = ObjectProperty(None)
    text_input = ObjectProperty(None)
    cancel = ObjectProperty(None)

def ScreenManager(ScreenManager):
    pass

class MenuScreen(Screen):
    
    container=None
    
    def download_ngs_primerplex(self):
        try:
            client=docker.from_env()
        except:
            print('ERROR (7)! Could not connect to docker VM')
            exit(7)
        t=threading.Thread(target=self.start_pulling,
                       args=(client,'aakechin/ngs-primerplex:native',))
        t.start()

    def prepare_reference(self):
        global global_image_names,available_images
        try:
            client=docker.from_env()
        except:
            print('ERROR (9)! Could not connect to docker VM')
            exit(9)
        t=threading.Thread(target=self.start_preparing_ref,
                           args=(client,global_image_names[1],))
        t.start()

    def start_pulling(self,client,image_name):
        client.images.pull(image_name)

    def start_preparing_ref(self,client,image_name):
        global global_image_names,available_images
        cmd=['gunzip',
             '/NGS-PrimerPlex/hg19/ucsc.hg19.fasta.gz',
             '/NGS-PrimerPlex/hg19/ucsc.hg19.fasta.amb.gz',
             '/NGS-PrimerPlex/hg19/ucsc.hg19.fasta.ann.gz',
             '/NGS-PrimerPlex/hg19/ucsc.hg19.fasta.bwt.gz',
             '/NGS-PrimerPlex/hg19/ucsc.hg19.fasta.fai.gz',
             '/NGS-PrimerPlex/hg19/ucsc.hg19.fasta.pac.gz',
             '/NGS-PrimerPlex/hg19/ucsc.hg19.fasta.sa.gz']
        imageFound=False
        try:
            images=client.images.list()
        except:
            print('WARNING (16)! Could not get list of available docker images!')
        for image in images:
            if len(image.tags)>0:
                try:
                    if global_image_names[0]==image.tags[0]:
                        imageFound=True
                        available_images[0]=True
                    if global_image_names[1]==image.tags[0]:
                        available_images[1]=True
                    if global_image_names[2]==image.tags[0]:
                        available_images[2]=True
                except:
                    print('WARNING (14)! Could not get tags for the following image:')
                    print(image.tags)
        if not imageFound:
            self.ids['prepare_ref_status']='ERROR! Download aakechin/ngs-primerplex:latest first'
        try:
            self.container=client.containers.run(image_name,
                                                 ' '.join(cmd),
                                                 detach=True,
                                                 entrypoint='')
        except:
            print('ERROR (15)! Could not run image:')
            print('Image name: '+image_name)
            print('Command: '+' '.join(cmd))
            exit(15)
        for stat in self.container.stats():
            self.ids['prepare_ref_status'].text+='.'
            if self.ids['prepare_ref_status'].text.count('.')>30:
                self.ids['prepare_ref_status'].text='.'
            containerStats=json.loads(stat.decode('utf-8'))
            if len(containerStats['pids_stats'])==0:
                break
        self.container.commit('ngs_primerplex_ref')
        self.container.remove()
        available_images[2]=True
        self.ids['prepare_ref_button'].disabled=True
        self.ids['prepare_ref_status'].text=''

    def on_enter(self,**kwargs):
        Clock.schedule_once(self.check_button_status)
        Clock.schedule_once(self.check_reference_status)

    def check_button_status(self,dt):
        global global_image_names,available_images
        try:
            client=docker.from_env()
        except:
            print('ERROR (8)! Could not connect to docker VM')
            exit(8)
        imageFound=False
        try:
            images=client.images.list()
        except:
            print('WARNING (13)! Could not get list of available docker images!')
        for image in images:
            if len(image.tags)>0:
                try:
                    if global_image_names[0]==image.tags[0]:
                        imageFound=True
                        available_images[0]=True
                    if global_image_names[1]==image.tags[0]:
                        imageFound=True
                        available_images[1]=True
                    if global_image_names[2]==image.tags[0]:
                        imageFound=True
                        available_images[2]=True
                except:
                    print('WARNING (14)! Could not get tags for the following image:')
                    print(image.tags)
        if not imageFound:
            self.ids['download_button'].disabled=False
            self.ids['extract_regions_menu_button'].disabled=True
            self.ids['design_primers_menu_button'].disabled=True
            self.ids['add_adapters_menu_button'].disabled=True
        else:
            self.ids['download_button'].disabled=True

    def check_reference_status(self,dt):
        global global_image_names,available_images
        try:
            client=docker.from_env()
        except:
            print('ERROR (11)! Could not connect to docker VM')
            exit(11)
        imageFound=False
        try:
            images=client.images.list()
        except:
            print('WARNING (12)! Could not get list of available docker images!')
        for image in images:
            if len(image.tags)>0:
                try:
                    if global_image_names[2]==image.tags[0]:
                        imageFound=True
                        available_images[2]=True
                except:
                    print('WARNING (10)! Could not get tags for the following image:')
                    print(image.tags)
        if not imageFound:
            self.ids['prepare_ref_button'].disabled=False
        else:
            self.ids['prepare_ref_button'].disabled=True

class SettingsScreen(Screen):
    def save_settings(self, path, filename):
        rFile=open(os.path.join(path, filename),'w')
        values={}
        for value in self.ids.keys():
            if value=='autoadjust':
                rFile.write('\t'.join([value,str(self.ids[value].active)])+'\n')
            else:
                rFile.write('\t'.join([value,str(self.ids[value].text)])+'\n')
        rFile.close()
        self.dismiss_popup()

    def load_settings(self, path, filename):
        file=open(os.path.join(path, filename[0]))
        for string in file:
            cols=string.replace('\n','').replace('\r','').split('\t')
            if cols[0]=='autoadjust' and cols[1]=='False':
                self.ids[cols[0]].active=False
            elif cols[0]=='autoadjust' and cols[1]=='True':
                self.ids[cols[0]].active=True
            else:
                self.ids[cols[0]].text=cols[1]
        file.close()
        self.dismiss_popup()
            
    def dismiss_popup(self):
        self._popup.dismiss()

    def show_load(self,title='Load file with settings'):
        content = LoadDialog(load=self.load_settings, cancel=self.dismiss_popup)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()

    def show_save(self,title='Choose file to save settings'):
        content = SaveDialog(save=self.save_settings, cancel=self.dismiss_popup)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()

class ExtractRegions(Screen):

    container=None
    
    def add_gene(self,*args):
        self.ids['genes_table'].rows+=1
        if round(self.ids['genes_table'].size_hint_y,1)<0.9:
            self.ids['genes_table'].size_hint_y+=0.1
            self.ids['genes_table_empty_label'].size_hint_y-=0.1
        if self.ids['genes_table'].children[0].size[1]<30:
            for textInput in self.ids['genes_table'].children:
                textInput.font_size-=0.5
        self.ids['genes_table'].add_widget(TextInput(multiline=False,
                                                     valign='top',
                                                     padding=(5,2,0,0)))
        self.ids['genes_table'].add_widget(TextInput(multiline=False,
                                                     valign='top',
                                                     padding=(5,2,0,0)))
        self.ids['genes_table'].add_widget(TextInput(multiline=False,
                                                     valign='top',
                                                     padding=(5,2,0,0)))

    def run_extraction(self,*args):
        global global_image_names,available_images
        app=App.get_running_app()
        self.ids['log_text_input'].text='Starting extraction of genome regions...'
        app.root.ids['menu_screen'].ids['extract_regions_menu_button'].background_color=(0.15,0.35,0.54,1)
        rows=[]
        # Reference directory
        refDir=self.ids['genbank_directory'].text
        if refDir=='<Chosen directory with GenBank-files>':
            refDir='/NGS-PrimerPlex/hg19/'
        else:
            if not os.path.isdir(refDir):
                refDir=os.path.dirname(refDir)
            if refDir[-1]!=os.path.sep:
                refDir+=os.path.sep
        # Reference genome
        wgRefPath=self.ids['wgref_file'].text
        if wgRefPath=='<Chosen reference genome FASTA-file>':
            wgRefPath='/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'
        wgRefDir=os.path.dirname(wgRefPath)
        wgRefName=os.path.basename(wgRefPath)
        # Output file
        outputFilePath=self.ids['output_file'].text
        if outputFilePath=='<Chosen output-file>':
            self.ids['log_text_input'].text='ERROR! Choose file for output!'
            return(None)
        outputFileDir=os.path.dirname(outputFilePath)
        outputFileName=os.path.basename(outputFilePath)
        # Read input table
        for i,cell in enumerate(self.ids['genes_table'].children):
            if i>=len(self.ids['genes_table'].children)-3:
                break
            elif i%3==0:
                rows.append([])
            elif i%3==2 and (cell.text=='' or
                             cell.text==' '):
                rows=rows[:-1]
            if len(rows)>0:
                rows[-1]=[cell.text.strip()]+rows[-1]
        if len(rows)==0:
            self.ids['log_text_input'].text='ERROR! Fill table with target genes!'
            return(None)
        file=open(outputFilePath+'_input_genes.txt','w')
        for row in rows:
            file.write('\t'.join(row)+'\n')
        file.close()
        # Intron size
        intronSize=self.ids['intron_size_text_input'].text
        # Command
        cmd=['python3',
             '/NGS-PrimerPlex/getGeneRegions.py',
             '-glf','/output/'+outputFileName+'_input_genes.txt',
             '-rf','/output/'+outputFileName,
             '-intron',intronSize]
        # Mounts and adding info about reference
        mounts=[]
        if refDir!='/NGS-PrimerPlex/hg19/':
            drive=os.path.splitdrive(refDir)
            if drive[0]!='':
                refDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
            mounts.append(Mount(target='/ref',
                                source=refDir,
                                type='bind'))
            cmd.extend(['-ref','/ref'])
            for image_name,available_image in zip(global_image_names,
                                                  available_images):
                # If image without reference is available, use it
                # Elif image with zipped hg19 available, use it
                if available_image:
                    break
        else:
            cmd.extend(['-ref',refDir])
            # In this case we can use only image with unzipped hg19
            if available_images[2]:
                image_name=global_image_names[2]
            elif available_images[1]:
                self.ids['log_text_input'].text='ERROR! Prepare hg19 reference first!'
                return(None)
            else:
                self.ids['log_text_input'].text='ERROR! Download aakechin/ngs-primerplex:latest first'
                return(None)
        drive=os.path.splitdrive(outputFileDir)
        if drive[0]!='':
            outputFileDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
        mounts.append(Mount(target='/output',
                            source=outputFileDir,
                            type='bind'))
        if wgRefPath!='/NGS-PrimerPlex/hg19/ucsc.hg19.fasta':
            drive=os.path.splitdrive(wgRefDir)
            if drive[0]!='':
                wgRefDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
            mounts.append(Mount(target='/wgref',
                                source=wgRefDir,
                                type='bind'))
            cmd.extend(['-wgref','/wgref/'+wgRefName])
        else:
            cmd.extend(['-wgref',wgRefPath])
            
        if self.ids['noncoding_switch'].active:
            cmd.append('-noncoding')
        t=threading.Thread(target=self.run_docker,
                           args=(image_name,cmd,mounts,))
        t.start()
        app.root.ids['design_primers'].ids['regions_file'].text=outputFilePath

    def run_docker(self,image_name,cmd,mounts):
        try:
            client=docker.from_env()
        except:
            print('ERROR (3)! Could not connect to docker VM')
            exit(3)
        try:
            self.container=client.containers.run(image_name,
                                                 ' '.join(cmd),
                                                 mounts=mounts,
                                                 detach=True,
                                                 entrypoint='')
        except:
            print('ERROR (4)! Could not run image:')
            print('Image name: '+image_name)
            print('Command: '+' '.join(cmd))
            print('Mounts: '+', '.join(map(str,mounts)))
            exit(4)
        for stat in self.container.stats():
            containerStats=json.loads(stat.decode('utf-8'))
            newLog=self.container.logs().decode('utf-8')
            if newLog!='':
                writtenLog=self.ids['log_text_input'].text
                if newLog!=writtenLog:
                    logParts=newLog.split('\n')
                    newLogParts=[]
                    for n,part in enumerate(logParts):
                        if '%' in part:
                            if n<len(logParts)-1 and '%' not in logParts[n+1]:
                                newLogParts.append(part)
                            elif n==len(logParts)-1:
                                newLogParts.append(part)
                        else:
                            newLogParts.append(part)
                    self.ids['log_text_input'].text='\n'.join(newLogParts)
            if len(containerStats['pids_stats'])==0:
                break
        # Change background of this step button on the main screen
        app=App.get_running_app()
        if 'NGS-PrimerPlex finished!' in self.container.logs().decode('utf-8'):
            app.root.ids['menu_screen'].ids['extract_regions_menu_button'].background_color=(0.226,0.527,0.273,1)
        else:
            app.root.ids['menu_screen'].ids['extract_regions_menu_button'].background_color=(0.754,0.02,0.226,1)
        self.container.remove()
        
    def dismiss_popup(self):
        self._popup.dismiss()

    def show_load(self,title='Load file'):
        if title=='Choose directory with GenBank-files':
            content = LoadDialog(load=self.loadRefDir, cancel=self.dismiss_popup)
        elif title=='Choose FASTA-file with reference genome':
            content = LoadDialog(load=self.loadRefGenome, cancel=self.dismiss_popup)
        else:
            print('ERROR (5)! Unknown file to load with the following title:')
            print(title)
            exit(5)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()

    def show_save(self,title='Output file'):
        content = SaveDialog(save=self.save, cancel=self.dismiss_popup)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()

    def loadRefDir(self, path, filename):
        self.ids['genbank_directory'].text=path
        self.dismiss_popup()

    def loadRefGenome(self, path, filename):
        self.ids['wgref_file'].text=os.path.join(path, filename[0])
        app=App.get_running_app()
        app.root.ids['design_primers'].ids['wgref_file'].text=os.path.join(path, filename[0])
        self.dismiss_popup()

    def save(self, path, filename):
        self.ids['output_file'].text=os.path.join(path, filename)
        self.dismiss_popup()

class DesignPrimers(Screen):

    container=None

    def start_design(self,*args):
        global global_image_names,available_images
        app=App.get_running_app()
        self.ids['design_primers_log'].text='Starting primer design...'
        app.root.ids['menu_screen'].ids['design_primers_menu_button'].background_color=(0.15,0.35,0.54,1)
        # Regions file
        regionsFilePath=self.ids['regions_file'].text
        if regionsFilePath=='<Choose file with target genome regions>':
            self.ids['design_primers_log'].text='ERROR! Choose file with target genome regions!'
            return(None)
        regionsFileDir=os.path.dirname(regionsFilePath)
        regionsFileName=os.path.basename(regionsFilePath)
        # Draft primers
        draftFilePath=self.ids['draft_file'].text
        if draftFilePath=='<Choose file with draft list of primers>':
            draftFilePath=None
        else:
            draftFileDir=os.path.dirname(draftFilePath)
            draftFileName=os.path.basename(draftFilePath)
        # Internal primers
        primersFilePath=self.ids['primers_file'].text
        if primersFilePath=='<Choose file with designed internal primers>':
            primersFilePath=None
        else:
            primersFileDir=os.path.dirname(primersFilePath)
            primersFileName=os.path.basename(primersFilePath)
        # Reference genome
        wgRefPath=self.ids['wgref_file'].text
        if wgRefPath=='<Chosen reference genome FASTA-file or leave it as it is to use hg19>':
            wgRefPath='/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'
        # VCF-file with SNPs
        snpsFilePath=self.ids['snps_file'].text
        if snpsFilePath=='<To check for covering SNPs choose VCF-file with SNPs>':
            snpsFilePath=None
        else:
            snpsFileDir=os.path.dirname(snpsFilePath)
            snpsFileName=os.path.basename(snpsFilePath)
        wgRefDir=os.path.dirname(wgRefPath)
        wgRefName=os.path.basename(wgRefPath)
        # Embedded PCR
        embedded=self.ids['embedded_switch'].active
        # Skip uncovered
        skip=self.ids['skip_switch'].active
        # Non-target hybridization
        nontargets=self.ids['nontargets_switch'].active
        # Embedded PCR
        snps=self.ids['snps_switch'].active
        # Threads number
        threads=self.ids['threads_text_input'].text.strip()
        # Run name
        runName=self.ids['run_name_text_input'].text.strip()
        drive=os.path.splitdrive(regionsFileDir)
        if drive[0]!='':
            regionsFileDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
        mounts=[Mount(target='/regions',
                      source=regionsFileDir,
                      type='bind')]
        cmd=['python3',
             '/NGS-PrimerPlex/NGS_primerplex.py',
             '-regions','/regions/'+regionsFileName,
             '-th',threads,
             '-run',runName]
        if wgRefPath!='/NGS-PrimerPlex/hg19/ucsc.hg19.fasta':
            drive=os.path.splitdrive(wgRefDir)
            if drive[0]!='':
                wgRefDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
            mounts.append(Mount(target='/wgref',
                                source=wgRefDir,
                                type='bind'))
            cmd.extend(['-ref','/wgref/'+wgRefName])
            for image_name,available_image in zip(global_image_names,
                                                  available_images):
                # If image without reference is available, use it
                # Elif image with zipped hg19 available, use it
                if available_image:
                    break
        else:
            cmd.extend(['-ref',wgRefPath])
            # In this case we can use only image with unzipped hg19
            if available_images[2]:
                image_name=global_image_names[2]
            elif available_images[1]:
                self.ids['design_primers_log'].text='ERROR! Prepare hg19 reference first!'
                return(None)
            else:
                self.ids['design_primers_log'].text='ERROR! Download aakechin/ngs-primerplex:latest first'
                return(None)
        if draftFilePath!=None:
            drive=os.path.splitdrive(draftFileDir)
            if drive[0]!='':
                draftFileDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
            mounts.append(Mount(target='/draft',
                                source=draftFileDir,
                                type='bind'))
            cmd.extend(['-draft','/draft/'+draftFileName])
        elif primersFilePath!=None:
            drive=os.path.splitdrive(primersFileDir)
            if drive[0]!='':
                primersFileDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
            mounts.append(Mount(target='/primers',
                                source=primersFileDir,
                                type='bind'))
            cmd.extend(['-primers','/primers/'+primersFileName])
        if snpsFilePath!=None:
            drive=os.path.splitdrive(snpsFileDir)
            if drive[0]!='':
                snpsFileDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
            mounts.append(Mount(target='/snps',
                                source=snpsFileDir,
                                type='bind'))
            cmd.extend(['-dbsnp','/snps/'+snpsFileName])
        if embedded:
            cmd.append('-embedded')
        if skip:
            cmd.append('-skip')
        if nontargets:
            cmd.append('-blast')
        if snps:
            cmd.append('-snps')
        cmd.append('-gui')
        if self.ids['left_adapter'].text.strip()!='':
            p=re.compile('([^ATGCatgc]+)')
            if len(p.findall(self.ids['left_adapter'].text.strip()))>0:
                self.ids['design_primers_log'].text='ERROR! Incorrect sequence of the left adapter!'
                return(None)
            else:
                cmd.extend(['-ad1',self.ids['left_adapter'].text.strip()])
        if self.ids['right_adapter'].text.strip()!='':
            p=re.compile('([^ATGCatgc]+)')
            if len(p.findall(self.ids['right_adapter'].text.strip()))>0:
                self.ids['design_primers_log'].text='ERROR! Incorrect sequence of the right adapter!'
                return(None)
            else:
                cmd.extend(['-ad2',self.ids['right_adapter'].text.strip()])
        # Values from settings
        for value in app.root.ids['settings'].ids.keys():
            if value!='autoadjust':
                cmd.extend(['-'+value,app.root.ids['settings'].ids[value].text])
            elif app.root.ids['settings'].ids[value].active==True:
                cmd.append('-'+value)
        self.ids['primer_design_stop_button'].disabled=False
        self.ids['primer_design_start_button'].disabled=True
        t=threading.Thread(target=self.run_docker,
                           args=(image_name,cmd,mounts,))
        t.start()
        # Change file name in the next step
        outputFileName=''.join([regionsFilePath[:-4],
                                '_NGS_primerplex_',
                                runName,
                                '_primers_combination_1_info.xls'])
        app.root.ids['add_adapters'].ids['filename_text_input_primers'].text=outputFileName

    def run_docker(self,image_name,cmd,mounts):
        try:
            client=docker.from_env()
        except:
            print('ERROR (3)! Could not connect to docker VM')
            exit(3)
        try:
            self.container=client.containers.run(image_name,
                                                 ' '.join(cmd),
                                                 mounts=mounts,
                                                 detach=True,
                                                 entrypoint='')
        except:
            print('ERROR (4)! Could not run image:')
            print('Image name: '+image_name)
            print('Command: '+' '.join(cmd))
            print('Mounts: '+', '.join(map(str,mounts)))
            exit(4)
        for stat in self.container.stats():
            containerStats=json.loads(stat.decode('utf-8'))
            newLog=self.container.logs().decode('utf-8')
            if newLog!='':
                writtenLog=self.ids['design_primers_log'].text
                if newLog!=writtenLog:
                    logParts=newLog.split('\n')
                    newLogParts=[]
                    for n,part in enumerate(logParts):
                        if '%' in part:
                            if n<len(logParts)-1 and '%' not in logParts[n+1]:
                                newLogParts.append(part)
                            elif n==len(logParts)-1:
                                newLogParts.append(part)
                        elif part=='':
                            continue
                        else:
                            newLogParts.append(part)
                    self.ids['design_primers_log'].text='\n'.join(newLogParts)
            if len(containerStats['pids_stats'])==0:
                break
        # Change background of this step button on the main screen
        app=App.get_running_app()
        if 'NGS-PrimerPlex finished!' in self.container.logs().decode('utf-8'):
            app.root.ids['menu_screen'].ids['design_primers_menu_button'].background_color=(0.226,0.527,0.273,1)
        else:
            app.root.ids['menu_screen'].ids['design_primers_menu_button'].background_color=(0.754,0.02,0.226,1)
        self.ids['primer_design_stop_button'].disabled=True
        self.ids['primer_design_start_button'].disabled=False
        self.container.remove()

    def stop_design(self):
        if self.container!=None:
            t=threading.Thread(target=self.stop_docker)
            t.start()            
        self.ids['primer_design_stop_button'].disabled=True
        self.ids['primer_design_start_button'].disabled=False
        self.ids['design_primers_log'].text+='\nSTOPPED'

    def stop_docker(self):
        self.container.stop()

    def dismiss_popup(self):
        self._popup.dismiss()

    def show_load(self,title='Load file'):
        app=App.get_running_app()
        if title=='Choose file with target genome regions':
            content = LoadDialog(load=self.loadTargetRegions, cancel=self.dismiss_popup)
        elif title=='Choose file with draft list of primers':
            content = LoadDialog(load=self.loadDraftPrimers, cancel=self.dismiss_popup)
        elif title=='Choose file with designed internal primers':
            content = LoadDialog(load=self.loadIntPrimers, cancel=self.dismiss_popup)
        elif title=='Choose FASTA-file with reference genome':
            content = LoadDialog(load=self.loadRefGenome, cancel=self.dismiss_popup)
        elif title=='Choose VCF-file with SNPs':
            content = LoadDialog(load=self.loadSnpsFile, cancel=self.dismiss_popup)
        else:
            print('ERROR (6)! Unknown file to load with the following title:')
            print(title)
            exit(6)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()

    def loadTargetRegions(self, path, filename):
        self.ids['regions_file'].text=os.path.join(path, filename[0])
        self.dismiss_popup()

    def loadDraftPrimers(self, path, filename):
        self.ids['draft_file'].text=os.path.join(path, filename[0])
        self.ids['primers_file'].text='<Choose file with designed internal primers>'
        self.dismiss_popup()

    def loadIntPrimers(self, path, filename):
        self.ids['primers_file'].text=os.path.join(path, filename[0])
        self.ids['draft_file'].text='<Choose file with draft list of primers>'
        self.ids['embedded_switch'].active=True
        self.dismiss_popup()

    def loadRefGenome(self, path, filename):
        self.ids['wgref_file'].text=os.path.join(path, filename[0])
        self.dismiss_popup()

    def loadSnpsFile(self, path, filename):
        self.ids['snps_file'].text=os.path.join(path, filename[0])
        self.ids['snps_switch'].active=True
        self.dismiss_popup()

class AddAdapters(Screen):

    container=None
    
    def dismiss_popup(self):
        self._popup.dismiss()

    def show_load_primers(self,title='Load file'):
        content = LoadDialog(load=self.load_primers, cancel=self.dismiss_popup)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()
        
    def load_primers(self, path, filename,):
        self.ids['filename_text_input_primers'].text=os.path.join(path,filename[0])
        self.dismiss_popup()

    def show_load_adapters(self,title='Load file'):
        content = LoadDialog(load=self.load_adapters, cancel=self.dismiss_popup)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()
        
    def load_adapters(self, path, filename,):
        self.ids['filename_text_input_adapters'].text=os.path.join(path,filename[0])
        self.dismiss_popup()
        
    def startAdding(self,*args):
        app=App.get_running_app()
        app.root.ids['menu_screen'].ids['add_adapters_menu_button'].background_color=(0.15,0.35,0.54,1)
        self.ids['add_adapters_log'].text='Started adding adapters...'
        t=threading.Thread(target=self.runAdding)
        t.start()
        
    def runAdding(self,*args):
        global global_image_names,available_images
        try:
            client=docker.from_env()
        except:
            print('ERROR (1)! Could not connect to docker VM')
            exit(1)
        inputFilePath=self.ids['filename_text_input_primers'].text
        inputDir=os.path.dirname(inputFilePath)
        inputFileName=os.path.basename(inputFilePath)
        drive=os.path.splitdrive(inputDir)
        if drive[0]!='':
            inputDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
        mounts=[Mount(target='/input',
                      source=inputDir,
                      type='bind')]
        cmd=['python3',
             '/NGS-PrimerPlex/addSeqToPrimers.py',
             '-in','/input/'+inputFileName]
        if self.ids['filename_text_input_adapters'].hint_text!='<Chosen file with adapters>':
            inputFilePath2=self.ids['filename_text_input_adapters'].text
            inputDir2=os.path.dirname(inputFilePath2)
            inputFileName2=os.path.basename(inputFilePath2)
            drive=os.path.splitdrive(inputDir2)
            if drive[0]!='':
                inputDir2='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
            mounts.append(Mount(target='/tags',
                                source=inputDir2,
                                type='bind'))
            cmd.extend(['-tags','/tags/'+inputFileName2])
        for image_name,available_image in zip(global_image_names,
                                              available_images):
            # If image without reference is available, use it
            # Elif image with zipped hg19 available, use it
            if available_image:
                break
        try:
            self.container=client.containers.run(image_name,
                                                 ' '.join(cmd),
                                                 mounts=mounts,
                                                 detach=True,
                                                 entrypoint='')
        except:
            print('ERROR (2)! Could not run image:')
            print('Image name: aakechin/ngs-primerplex')
            print('Command: '+' '.join(cmd))
            print('Mounts: '+', '.join(map(str,mounts)))
            exit(2)
        for stat in self.container.stats():
            containerStats=json.loads(stat.decode('utf-8'))
            newLog=self.container.logs().decode('utf-8')
            if newLog!='':
                writtenLog=self.ids['add_adapters_log'].text
                if newLog!=writtenLog:
                    logParts=newLog.split('\n')
                    newLogParts=[]
                    for n,part in enumerate(logParts):
                        if '%' in part:
                            if n<len(logParts)-1 and '%' not in logParts[n+1]:
                                newLogParts.append(part)
                            elif n==len(logParts)-1:
                                newLogParts.append(part)
                        else:
                            newLogParts.append(part)
                    self.ids['add_adapters_log'].text='\n'.join(newLogParts)
            if len(containerStats['pids_stats'])==0:
                break
        app=App.get_running_app()
        app.root.ids['menu_screen'].ids['add_adapters_menu_button'].background_color=(0.226,0.527,0.273,1)
        self.container.remove()

class ConvertToDraft(Screen):

    container=None
    
    def dismiss_popup(self):
        self._popup.dismiss()

    def show_load_primers(self,title='Load file'):
        content = LoadDialog(load=self.load_primers, cancel=self.dismiss_popup)
        self._popup = Popup(title=title, content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()
        
    def load_primers(self, path, filename,):
        self.ids['filename_text_input_primers'].text=os.path.join(path,filename[0])
        self.dismiss_popup()
        
    def startConverting(self,*args):
        app=App.get_running_app()
        app.root.ids['menu_screen'].ids['convert_to_draft_menu_button'].background_color=(0.15,0.35,0.54,1)
        self.ids['convert_to_draft_log'].text='Started converting to draft-primers...'
        t=threading.Thread(target=self.runAdding)
        t.start()
        
    def runAdding(self,*args):
        global global_image_names,available_images
        try:
            client=docker.from_env()
        except:
            print('ERROR (1)! Could not connect to docker VM')
            exit(1)
        inputFilePath=self.ids['filename_text_input_primers'].text
        inputDir=os.path.dirname(inputFilePath)
        inputFileName=os.path.basename(inputFilePath)
        drive=os.path.splitdrive(inputDir)
        if drive[0]!='':
            inputDir='/'+drive[0].replace(':','')+drive[1].replace('\\','/')
        mounts=[Mount(target='/input',
                      source=inputDir,
                      type='bind')]
        cmd=['python3',
             '/NGS-PrimerPlex/convertToDraftFile.py',
             '-in','/input/'+inputFileName,
             '-out','/input/'+inputFileName[:inputFileName.rfind('.')]+'.draft.xls']
        for image_name,available_image in zip(global_image_names,
                                              available_images):
            # If image without reference is available, use it
            # Elif image with zipped hg19 available, use it
            if available_image:
                break
        try:
            self.container=client.containers.run(image_name,
                                                 ' '.join(cmd),
                                                 mounts=mounts,
                                                 detach=True,
                                                 entrypoint='')
        except:
            print('ERROR (2)! Could not run image:')
            print('Image name: aakechin/ngs-primerplex')
            print('Command: '+' '.join(cmd))
            print('Mounts: '+', '.join(map(str,mounts)))
            exit(2)
        for stat in self.container.stats():
            containerStats=json.loads(stat.decode('utf-8'))
            newLog=self.container.logs().decode('utf-8')
            if newLog!='':
                writtenLog=self.ids['convert_to_draft_log'].text
                if newLog!=writtenLog:
                    logParts=newLog.split('\n')
                    newLogParts=[]
                    for n,part in enumerate(logParts):
                        if '%' in part:
                            if n<len(logParts)-1 and '%' not in logParts[n+1]:
                                newLogParts.append(part)
                            elif n==len(logParts)-1:
                                newLogParts.append(part)
                        else:
                            newLogParts.append(part)
                    self.ids['convert_to_draft_log'].text='\n'.join(newLogParts)
            if len(containerStats['pids_stats'])==0:
                break
        app=App.get_running_app()
        app.root.ids['menu_screen'].ids['convert_to_draft_menu_button'].background_color=(0.226,0.527,0.273,1)
        self.container.remove()

class NGS_PrimerPlexApp(App):
    pass

if __name__ == '__main__':
    NGS_PrimerPlexApp().run()
