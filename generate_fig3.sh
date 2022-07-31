ts=(1.45 36)

for THRESHOLD in ${ts[@]};
do
    ### Inputs
    numberOfProcessors=2
    numberOfSimulations=2

    coralPercent=33
    macroalgaePercent=33
    gridOption=0 #array

    #Time and Grid Settings
    rows=25
    columns=25
    neighborhoodThreshold=$THRESHOLD

    recordRate=10
    imageReturn=true #lowercase spelling true or false
    imageRecordRate=1 ## after how many recordings do you want an image saved? 1 = each time the other recordings are taken

        ## table for first run. later averages
    #loop through g, nested loop through gridoption (once ready)
    r=1.0
    d=.4
    a=.2
    g=.54 #array
    y=.75
    dt=.1
    tf=11 #can play with this value as well

    blobValue=0

    ## create directory
   grazingFolder=$(python -c "print(str(round($g*100)))")
   thresholdFolder=$(python -c "print(int($neighborhoodThreshold*100))")


    mkdir -p 'output'/$rows'x'$columns/'grid'$gridOption/'grazing'$grazingFolder/'threshold'$thresholdFolder

    ## run simulation
    python coralModel_functions.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $imageRecordRate $r $d $a $g $y $dt $tf $blobValue

done
