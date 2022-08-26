# First compute landscapes for all grazing rates
numberOfProcessors=100
numberOfSimulations=100
coralPercent=33
macroalgaePercent=33
rows=25
columns=25
neighborhoodThreshold=1.45
gridOption=0
recordRate=10
imageReturn=0
imageRecordRate=1 ## after how many recordings do you want an image saved? 1 = each time the other recordings are taken
r=1.0
d=0.4
a=.2
#g
y=.75
dt=.1
blobValue=0
tf=101
organism=0
group=1
denoise1=1
denoise2=99
isUnion=0
isCoral=1

# g values to use
gs=(0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63)

for g in ${gs[@]};
do
                        

    grazingFolder=$(python -c "print((round(100*$g)))")
    thresholdFolder=$(python -c "print(int($neighborhoodThreshold*100))")
    
    tims=$tf
    # run simulation
        
    # create directory for bars
    mkdir -p 'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/

    # generate bars
    python PL_grazing_print_barcodes.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $imageRecordRate $r $d $a $g $y $dt $tf $organism $group $denoise1 $denoise2 $isUnion $tims

    # convert to landscapes
    xargs <./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'c'$coralPercent'm'$macroalgaePercent'all_bars_list_hdim0.txt' -I filename ./gudhi.3.3.0/build/utilities/Persistence_representations/persistence_landscapes/create_landscapes -1 "filename"

    # compute average landscape
    ./gudhi.3.3.0/build/utilities/Persistence_representations/persistence_landscapes/average_landscapes $(cat ./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'c'$coralPercent'm'$macroalgaePercent'all_lands_list_hdim0.txt')

    #mv average.land ./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/

done

# Now run simulations in the metastable region for longer

isUnion=0
isCoral=1
isMacro=0

gs=(0.52 0.53 0.54 0.55)

for g in ${gs[@]};
do

    ## create directory
    grazingFolder=$(python -c "print(round($g*100))")
    thresholdFolder=$(python -c "print(round($neighborhoodThreshold*100))")

    tf=1000
    tims=101

 
    ## generate barcodes
    mkdir -p 'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'coral_success'
    mkdir -p 'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'macro_success'

    python PL_prediction_print_barcodes.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $recordRate $r $d $a $g $y $dt $tf $organism $group $denoise1 $denoise2 $isUnion $tims

    # For coral success

    # convert to landscapes
    xargs <./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'coral_success'/'c'$coralPercent'm'$macroalgaePercent'coral_success_all_hdim0.txt' -I filename ./gudhi.3.3.0/build/utilities/Persistence_representations/persistence_landscapes/create_landscapes -1 "filename"

    # compute average landscape
    ./gudhi.3.3.0/build/utilities/Persistence_representations/persistence_landscapes/average_landscapes $(cat ./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'coral_success'/'c'$coralPercent'm'$macroalgaePercent'coral_success_all_land_hdim0.txt')
    mv average.land ./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'coral_success'/


    # For macro success

    # convert to landscapes
    xargs <./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'macro_success'/'c'$coralPercent'm'$macroalgaePercent'macro_success_all_hdim0.txt' -I filename ./gudhi.3.3.0/build/utilities/Persistence_representations/persistence_landscapes/create_landscapes -1 "filename"

    # compute average landscape
    ./gudhi.3.3.0/build/utilities/Persistence_representations/persistence_landscapes/average_landscapes $(cat ./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'macro_success'/'c'$coralPercent'm'$macroalgaePercent'macro_success_all_land_hdim0.txt')
    mv average.land ./'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder/'bars'/'macro_success'/

done
