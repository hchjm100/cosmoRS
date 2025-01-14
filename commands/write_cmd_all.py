import os


file_name = 'run2'
snapshot_number = 62
folder_name = os.listdir('../../') 
folder_name = 'Gadget2/SIDM/codeSIDM3/output_testRun2_job'
#     [folder for folder in folder_name if 'DM' in folder]
# folder_name.sort()
print(folder_name)

with open(file_name + '.sh', 'w') as sh:
    sh.write('rm ../../' + folder_name + '/boundmass.txt\n')
    sh.write('rm ../../' + folder_name + '/dataProfile/*\n')
    for i in range(snapshot_number):
        sh.write('nohup srun --mem=10gb --pty ./modRS ../../' + folder_name + ' snapshot {:03} -l &\n'.format(i))
sh.close()

# j=0

# with open("sh-all.sh", 'w') as runall:
#   for folder in folder_name:
#       runall.write('sh ' + folder + '.sh\n')
#       j+=1
#       if j%2==0:
#         runall.write('sleep 135\n')
#       with open(folder + '.sh','w') as sh:
#          sh.write('rm ../../'+folder+'/output/boundmass.txt\n')
#          sh.write('rm ../../'+folder+'/output/dataProfile/*\n')
#          for i in range(60):
#             sh.write('nohup srun --mem=10gb --pty ./modRS ../../' +folder+ '/output snapshot {:03} -l &\n'.format(i))
#       sh.close()
# runall.close




