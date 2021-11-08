read -p "Please tell me the protein name you are interested in: " pro
 
read -p "Do you know which taxonomic group it is?(y/n) " a
if test $a == 'y';then 
    read -p "Pleas tell me the name of its taxonomic group: " group
    else
    $group=True
fi

db="protein"

# Set the query key and WebEnv
wget -qO- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=$db&term=$beast&usehistory=y" | awk '{
if(FNR==3)
    {
        split($0,xml_line_split,"<|>");
        query_key=xml_line_split[17];
        WebEnv=xml_line_split[21];
        print query_key > "query_key";
        print WebEnv > "WebEnv";
        exit;}}'

        query_key=$(cat query_key)
        WebEnv=$(cat WebEnv)


# Use esearch then efetch to get the UID values
echo "Finding the uids that related to $pro"
if test $group ;
    then 
    esearch -db $db -query $pro | efetch -format uid > aim.uid
    else 
    esearch -db $db -query $pro AND $group | efetch -format uid > aim.uid
fi
