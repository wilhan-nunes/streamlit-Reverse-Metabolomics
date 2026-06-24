import json
import time
import logging

import pandas as pd
import argparse
import requests
import requests_cache
from tqdm import tqdm

HOST = "https://api.fasst.gnps2.org"


def _blocking_for_results(task_id: str, host: str = HOST, retries_max: int = 120) -> dict:
    for _ in range(retries_max):
        r = requests.get(f"{host}/search/result/{task_id}", timeout=30)
        r.raise_for_status()
        payload = r.json()
        if isinstance(payload, dict) and payload.get("status") == "PENDING":
            time.sleep(1)
            continue
        return payload
    raise TimeoutError("Timeout waiting for results from FASST API")


def masst_query_all(query_df, database, analog=False, lower_delta=130, upper_delta=200, precursor_mz_tol=0.02, fragment_mz_tol=0.02, min_cos=0.7):
    output_results_list = []

    for query_element in tqdm(query_df.to_dict(orient="records")):
        try:
            print(query_element)
            usi = query_element["usi"]

            params = {
                "library": str(database),
                "usi": usi,
                "analog": "Yes" if analog else "No",
                "lower_delta": lower_delta,
                "upper_delta": upper_delta,
                "pm_tolerance": precursor_mz_tol,
                "fragment_tolerance": fragment_mz_tol,
                "cosine_threshold": min_cos,
            }

            search_api_response = requests.post(f"{HOST}/search", json=params, timeout=30)
            search_api_response.raise_for_status()
            task_id = search_api_response.json()["id"]

            search_api_response_json = _blocking_for_results(task_id)

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
    