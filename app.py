import streamlit as st
import pandas as pd
from rdkit import Chem
from io import BytesIO

st.set_page_config(page_title="Filtro SMARTS", layout="wide")

st.title("🔬 Filtro de moléculas por SMARTS")

# =========================
# INGRESO SMARTS
# =========================

user_smarts = st.text_input(
    "Ingrese SMARTS (puede colocar varios separados por coma o salto de línea)",
    placeholder="Ejemplo: c1ccc2ncnc2c1"
)

# =========================
# SELECCIONAR DELIMITADOR
# =========================

st.sidebar.header("Configuración del archivo")

separator_option = st.sidebar.radio(
    "Seleccione el separador del CSV:",
    options=[";", ","],
    index=0  # Por defecto punto y coma
)

st.sidebar.write(f"Separador seleccionado: '{separator_option}'")

# =========================
# CARGA ARCHIVO
# =========================

archivo = st.file_uploader(
    "Sube tu archivo CSV (debe contener la columna 'SMILES')",
    type=["csv"]
)

if archivo is not None:

    # Leer usando el separador seleccionado
    df = pd.read_csv(archivo, sep=separator_option)

    # Limpiar posibles espacios
    df.columns = df.columns.str.strip()

    if "SMILES" not in df.columns:
        st.error(f"La columna 'SMILES' no fue encontrada. Columnas detectadas: {list(df.columns)}")
        st.stop()

    st.success("Archivo cargado correctamente")

    if user_smarts.strip() == "":
        st.warning("Ingrese al menos un SMARTS")
        st.stop()

    # =========================
    # PROCESAR SMARTS
    # =========================

    smarts_list = []
    raw_smarts = user_smarts.replace(",", "\n").split("\n")

    for smarts in raw_smarts:
        smarts = smarts.strip()
        if smarts != "":
            smarts_list.append(smarts)

    patterns = []
    for smarts in smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            st.error(f"SMARTS inválido: {smarts}")
            st.stop()
        patterns.append(patt)

    # =========================
    # FUNCION MATCH
    # =========================

    def match_any_smarts(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        for patt in patterns:
            if mol.HasSubstructMatch(patt):
                return True
        return False

    # =========================
    # EJECUTAR
    # =========================

    if st.button("🔎 Ejecutar búsqueda"):

        with st.spinner("Analizando moléculas..."):

            df["match"] = df["SMILES"].apply(match_any_smarts)

            df_matches = df[df["match"]].drop(columns=["match"])
            df_filtered = df[~df["match"]].drop(columns=["match"])

            indices_matches = df[df["match"]].index.tolist()

        # =========================
        # RESULTADOS
        # =========================

        st.subheader("📊 Resultados")

        col1, col2, col3 = st.columns(3)

        col1.metric("Moléculas totales", len(df))
        col2.metric("Coincidencias", len(df_matches))
        col3.metric("Después del filtrado", len(df_filtered))

        # =========================
        # INDICES
        # =========================

        st.subheader("📌 Índices coincidentes")
        st.code(str(indices_matches), language="python")

        # =========================
        # DESCARGAS
        # =========================

        def to_csv_bytes(dataframe):
            buffer = BytesIO()
            dataframe.to_csv(buffer, index=True, sep=separator_option)
            return buffer.getvalue()

        def indices_to_csv(indices):
            df_idx = pd.DataFrame({"index": indices})
            buffer = BytesIO()
            df_idx.to_csv(buffer, index=False, sep=separator_option)
            return buffer.getvalue()

        st.subheader("⬇️ Descargar")

        st.download_button(
            "(a) Descargar SOLO coincidencias",
            data=to_csv_bytes(df_matches),
            file_name="coincidencias_smarts.csv",
            mime="text/csv"
        )

        st.download_button(
            "(b) Descargar base SIN coincidencias",
            data=to_csv_bytes(df_filtered),
            file_name="base_filtrada.csv",
            mime="text/csv"
        )

        st.download_button(
            "Descargar índices coincidentes",
            data=indices_to_csv(indices_matches),
            file_name="indices_coincidentes.csv",
            mime="text/csv"
        )
