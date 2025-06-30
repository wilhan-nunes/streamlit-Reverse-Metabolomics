import pandas as pd
import argparse
import requests
import requests_cache
from tqdm import tqdm
#requests_cache.install_cache('demo_cache')

def query_usi(usi, database, analog=False, precursor_mz_tol=0.02, fragment_mz_tol=0.02, min_cos=0.7):
    URL = "https://fasst.gnps2.org/search"

    params = {
        "usi": usi,
        "library": database,
        "analog": "Yes" if analog else "No",
        "pm_tolerance": precursor_mz_tol,
        "fragment_tolerance": fragment_mz_tol,
        "cosine_threshold": min_cos,
    }

    r = requests.get(URL, params=params, timeout=50)
    
    r.raise_for_status()

    return r.json()

def masst_query_all(query_df, database, masst_type, analog=False, precursor_mz_tol=0.02, fragment_mz_tol=0.02, min_cos=0.7):
    output_results_list = []

    for query_element in tqdm(query_df.to_dict(orient="records")):
        try:

            usi = query_element["usi"]

            results_dict = query_usi(usi, database,
                analog=analog, precursor_mz_tol=precursor_mz_tol, 
                fragment_mz_tol=fragment_mz_tol, min_cos=min_cos)
            results_df = pd.DataFrame(results_dict["results"])

            # TODO: Support munging of microbemasst results
            #if masst_type == "microbemasst":
                # Lets do additionally processing
            #    print("MICROBEMASST")
            # TODO: Merge with metadata automatically

            results_df["query_usi"] = usi
            if "flag" in query_element:
                results_df["flag"] = query_element["flag"]

            output_results_list.append(results_df)
        except:
            pass
    
    output_results_df = pd.concat(output_results_list)

    return output_results_df

def main():
    parser = argparse.ArgumentParser(description='Fast MASST Client')
    parser.add_argument('input_file', help='file to query with USIs')
    parser.add_argument('output_file', help='output_file')
    parser.add_argument('--masst_type', help='Type of MASST to give youresults: gnpsdata, microbemasst', default="masst")
    parser.add_argument('--database', help='Type database to actually search', default="gnpsdata_index")
    parser.add_argument('--analog', help='Perform Yes or No', default="No")
    parser.add_argument('--precursor_tolerance', help='precursor_tolerance', default=0.02, type=float)
    parser.add_argument('--fragment_tolerance', help='fragment_tolerance', default=0.02, type=float)
    parser.add_argument('--cosine', help='cosine', default=0.7, type=float)

    args = parser.parse_args()

    analog_boolean = args.analog == "Yes"

    query_df = pd.read_csv(args.input_file, sep=None)

    output_results_df = masst_query_all(query_df, 
                                        args.database, args.masst_type, 
                                        analog=analog_boolean,
                                        precursor_mz_tol=args.precursor_tolerance,
                                        fragment_mz_tol=args.fragment_tolerance)
                                        
    output_results_df.to_csv(args.output_file, index=False, sep="\t")

if __name__ == '__main__':
    main()
    