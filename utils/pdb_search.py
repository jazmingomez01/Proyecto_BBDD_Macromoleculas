import re
import requests
import time
import pandas as pd

from data.fetch_uniprot import buscar_pdb_uniprot


# Verifica si el accession es de UniProt
def es_accession_uniprot(accession):
    return bool(
        re.match(
            r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[A-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[A-Z0-9]{7}$",
            accession,
        )
    )


# Formatea los resultados de PDB en una tabla usando pandas
# Entrada = resultados de PDB
# Salida = tabla formateada con información de las estructuras PDB
def formatear_resultados_pdb(resultados):

    if isinstance(resultados, str):
        return resultados  # Es un mensaje de error

    if not resultados:
        return "No se encontraron estructuras PDB"

    # Crear DataFrame con la nueva estructura
    df = pd.DataFrame(resultados)

    # Renombrar columnas para mejor presentación
    df = df.rename(
        columns={
            "identifier": "Identifier",
            "method": "Method",
            "resolution": "Resolution",
            "chain": "Chain/Position",
        }
    )

    # Configurar pandas para mejor visualización
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", 50)

    return df.to_string(index=False, justify="left")


# Busca estructuras PDB asociadas a un accession
# Entrada = accession de UniProt
# Salida = tabla con información detallada de PDBs encontrados
def lista_pdb(accession):

    print(f"Buscando estructuras PDB para: {accession}")

    #Resolver accession → lista de UniProt IDs
    if es_accession_uniprot(accession):
        uniprot_ids = [accession]
    else:
        try:
            uniprot_ids = map_ncbi_to_uni(accession)
            if not uniprot_ids:
                return f"Error: '{accession}' no es un accession válido de UniProt"
        except Exception as e:
            return f"Error durante el mapeo: {e}"

    #Buscar PDBs para cada UniProt ID
    total_pdbs = []
    for uid in uniprot_ids:
        pdb_info = buscar_pdb_uniprot(uid)
        if pdb_info:
            total_pdbs.extend(pdb_info)

    if not total_pdbs:
        return f"No se encontraron estructuras PDB para el ID '{accession}' (UniProt: {', '.join(uniprot_ids)})"

    # Usar pandas para formatear y mostrar los resultados
    tabla_formateada = formatear_resultados_pdb(pdb_info)

    print(f"\nEncontradas {len(pdb_info)} estructuras PDB:")
    print("=" * 120)
    print(tabla_formateada)
    print("=" * 120)

# Mapea un ID de NCBI a UniProt usando la API de UniProt
# Entrada = ID de NCBI (RefSeq Protein)
# Salida = lista de IDs de UniProt correspondientes
# Excepciones: Lanza RuntimeError si hay problemas de conexión o el mapeo
def map_ncbi_to_uni(ncbi_id):
    run_url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "RefSeq_Protein",
        "to": "UniProtKB",
        "ids": ncbi_id
    }

    try:
        response = requests.post(run_url, data=params, timeout=30)
        response.raise_for_status()
        job_id = response.json()["jobId"]
    except Exception as e:
        raise RuntimeError(f"Error al enviar la solicitud de mapeo: {e}")

    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"

    while True:
        try:
            time.sleep(3)
            status_response = requests.get(status_url, timeout=30)
            status_response.raise_for_status()
            status_data = status_response.json()

            if "jobStatus" in status_data:
                job_status = status_data["jobStatus"]
                if job_status == "FINISHED":
                    break
                elif job_status == "FAILED":
                    raise RuntimeError("El mapeo falló.")
            elif "results" in status_data:
                break
        except Exception as e:
            raise RuntimeError(f"Error al verificar el estado del trabajo: {e}")

    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"

    try:
        result_response = requests.get(result_url, timeout=30)
        result_response.raise_for_status()
        results_data = result_response.json()
    except Exception as e:
        raise RuntimeError(f"Error al obtener los resultados: {e}")

    uniprot_ids = [item["to"] for item in results_data.get("results", [])]

    return uniprot_ids
