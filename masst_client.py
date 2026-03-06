import json

import pandas as pd
import argparse
import requests
import requests_cache
from tqdm import tqdm
from gnpsdata import fasst
#requests_cache.install_cache('demo_cache')

URL = "https://fasst.gnps2.org/search"

def masst_query_all(query_df, database, analog=False, lower_delta=130, upper_delta=200, precursor_mz_tol=0.02, fragment_mz_tol=0.02, min_cos=0.7):
    output_results_list = []

    for query_element in tqdm(query_df.to_dict(orient="records")):
        try:
            print(query_element)
            usi = query_element["usi"]

            #retrieve usi json from metabolomics resolver API
            usi_json = requests.get(f"https://metabolomics-usi.gnps2.org/json/?usi1={usi}").json()
            # modify precursor_charge value to be an absolute value (since some usis have negative charge values)
            usi_json["precursor_charge"] = abs(usi_json["precursor_charge"])
            # spec_json = json.dumps(usi_json)
    
            dps = [[round(mz, 5), intensity] for mz, intensity in usi_json["peaks"]]
            spec_dict = {
                "n_peaks": len(dps),
                "peaks": dps,
                "precursor_mz": usi_json["precursor_mz"],
                "precursor_charge": abs(usi_json["precursor_charge"]),
            }

            max_intensity = max([v[1] for v in spec_dict["peaks"]])
            dps = [
                [round(dp[0], 5), round(dp[1] / max_intensity * 100.0, 1)]
                for dp in spec_dict["peaks"]
            ]
            dps = [dp for dp in dps if dp[1] >= 0.1]

            spec_dict["peaks"] = dps
            spec_dict["n_peaks"] = len(dps)
    
            spec_json = json.dumps(spec_dict)
    
            params = {
            "library": str(database),
            "analog": "Yes" if analog else "No",
            "delta_mass_below": lower_delta,
            "delta_mass_above": upper_delta,
            "pm_tolerance": precursor_mz_tol,
            "fragment_tolerance": fragment_mz_tol,
            "cosine_threshold": min_cos,
            "query_spectrum": spec_json,
        }

            search_api_response = requests.post(URL, data=params, timeout=300)
            search_api_response.raise_for_status()
            search_api_response_json = search_api_response.json()

            results_df = pd.DataFrame(search_api_response_json['results'])

            results_df["query_usi"] = usi
            if "flag" in query_element:
                results_df["flag"] = query_element["flag"]

            output_results_list.append(results_df)
        except:
            print(f"Error occurred while processing USI: {usi}")

    if len(output_results_list) == 0:
        return pd.DataFrame()
    output_results_df = pd.concat(output_results_list)

    return output_results_df

def main():
    parser = argparse.ArgumentParser(description='Fast MASST Client')
    parser.add_argument('input_file', help='file to query with USIs')
    parser.add_argument('output_file', help='output_file')
    parser.add_argument('--database', help='Type database to actually search', default="gnpsdata_index")
    parser.add_argument('--analog', help='Perform Yes or No', default="No")
    parser.add_argument('--lower_delta', help='Lower Precursor m/z Delta (Da) for analog search', default=130.0, type=float)
    parser.add_argument('--upper_delta', help='Upper Precursor m/z Delta (Da) for analog search', default=200.0, type=float)
    parser.add_argument('--precursor_tolerance', help='precursor_tolerance', default=0.02, type=float)
    parser.add_argument('--fragment_tolerance', help='fragment_tolerance', default=0.02, type=float)
    parser.add_argument('--cosine', help='cosine', default=0.7, type=float)

    args = parser.parse_args()

    analog_boolean = args.analog == "Yes"

    query_df = pd.read_csv(args.input_file, sep=None)

    output_results_df = masst_query_all(query_df, 
                                        args.database, 
                                        analog=analog_boolean,
                                        lower_delta=args.lower_delta,
                                        upper_delta=args.upper_delta,
                                        precursor_mz_tol=args.precursor_tolerance,
                                        fragment_mz_tol=args.fragment_tolerance)
                                        
    output_results_df.to_csv(args.output_file, index=False, sep="\t")

if __name__ == '__main__':
    main()
    