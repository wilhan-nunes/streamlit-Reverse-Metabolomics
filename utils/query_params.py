"""
Query parameter utilities for URL-based configuration sharing.

Schema-based approach for maintainability and external tool integration.
"""

import json
import urllib.parse
from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, TypeVar, Generic
import pandas as pd
import streamlit as st


T = TypeVar('T')


class ParamType(Enum):
    """Supported parameter types for external documentation."""
    STRING = "string"
    BOOLEAN = "boolean"
    FLOAT = "float"
    INTEGER = "integer"
    JSON = "json"  # For complex objects like USI dict


@dataclass
class QueryParam(Generic[T]):
    """
    Definition of a single query parameter.
    
    Attributes:
        key: URL parameter key
        default: Default value if not in URL
        param_type: Type for documentation/validation
        description: Human-readable description for external tools
        encoder: Function to encode value to URL string
        decoder: Function to decode URL string to value
        validator: Optional function to validate decoded value
        choices: Optional list of valid values (for enums/selects)
    """
    key: str
    default: T
    param_type: ParamType
    description: str = ""
    encoder: Callable[[T], str] = None
    decoder: Callable[[str], T] = None
    validator: Callable[[T], bool] = None
    choices: list[Any] = None
    
    def __post_init__(self):
        # Set default encoder/decoder based on type
        if self.encoder is None:
            self.encoder = self._default_encoder()
        if self.decoder is None:
            self.decoder = self._default_decoder()
    
    def _default_encoder(self) -> Callable[[T], str]:
        if self.param_type == ParamType.BOOLEAN:
            return lambda v: str(v).lower()
        elif self.param_type == ParamType.JSON:
            return lambda v: urllib.parse.quote(json.dumps(v, separators=(',', ':')), safe='')
        else:
            return str
    
    def _default_decoder(self) -> Callable[[str], T]:
        if self.param_type == ParamType.BOOLEAN:
            return lambda v: str(v).lower() in ('true', '1', 'yes')
        elif self.param_type == ParamType.FLOAT:
            return lambda v: float(v)
        elif self.param_type == ParamType.INTEGER:
            return lambda v: int(v)
        elif self.param_type == ParamType.JSON:
            return lambda v: json.loads(urllib.parse.unquote(v))
        else:
            return str
    
    def encode(self, value: T) -> str | None:
        """Encode value to URL string. Returns None if equals default."""
        # Handle DataFrame comparison specially (can't use == with None)
        if isinstance(value, pd.DataFrame):
            if value.empty:
                return None
        elif value is None and self.default is None:
            return None
        elif value == self.default:
            return None
        
        try:
            encoded = self.encoder(value)
            # Return None for empty strings
            if encoded == "" or encoded is None:
                return None
            return encoded
        except (ValueError, TypeError):
            return None
    
    def decode(self, raw: str) -> T:
        """Decode URL string to value. Returns default on failure."""
        try:
            value = self.decoder(raw)
            if self.validator and not self.validator(value):
                return self.default
            if self.choices and value not in self.choices:
                return self.default
            return value
        except (ValueError, TypeError, json.JSONDecodeError):
            return self.default
    
    def to_schema(self) -> dict:
        """Export parameter definition for external tools (OpenAPI-style)."""
        schema = {
            "name": self.key,
            "type": self.param_type.value,
            "default": self.default,
            "description": self.description,
        }
        if self.choices:
            schema["enum"] = [c for c in self.choices if c is not None]
        return schema


# =============================================================================
# PARAMETER REGISTRY - Add new parameters here
# =============================================================================

PARAM_REGISTRY: dict[str, QueryParam] = {}


def register_param(param: QueryParam) -> QueryParam:
    """Register a parameter in the global registry."""
    PARAM_REGISTRY[param.key] = param
    return param


# --- USI Data (special handling) ---
def encode_usi_df(usi_df: pd.DataFrame) -> str:
    """Encode USI DataFrame to URL-safe JSON."""
    if usi_df is None or len(usi_df) == 0:
        return ""
    valid_df = usi_df[
        (usi_df['usi'].str.strip() != '') & 
        (~usi_df['usi'].str.contains('^USI-\\d+$', regex=True, na=False)) &
        (usi_df['compound_name'].str.strip() != '') &
        (~usi_df['compound_name'].str.contains('^Compound-name-\\d+$', regex=True, na=False))
    ]
    if len(valid_df) == 0:
        return ""
    usi_dict = dict(zip(valid_df['compound_name'], valid_df['usi']))
    return urllib.parse.quote(json.dumps(usi_dict, separators=(',', ':')), safe='')


def decode_usi_df(encoded: str) -> pd.DataFrame | None:
    """Decode URL-encoded JSON to USI DataFrame."""
    if not encoded:
        return None
    try:
        usi_dict = json.loads(urllib.parse.unquote(encoded))
        if not isinstance(usi_dict, dict) or len(usi_dict) == 0:
            return None
        return pd.DataFrame({
            'usi': list(usi_dict.values()),
            'compound_name': list(usi_dict.keys())
        })[['usi', 'compound_name']]
    except (json.JSONDecodeError, ValueError, TypeError):
        return None


# --- Register all parameters ---

USI_PARAM = register_param(QueryParam(
    key="usi",
    default=None,
    param_type=ParamType.JSON,
    description="USI data as JSON object: {\"compound_name\": \"mzspec:...\", ...}",
    encoder=encode_usi_df,
    decoder=decode_usi_df,
))

EXAMPLE_PARAM = register_param(QueryParam(
    key="example",
    default=None,
    param_type=ParamType.STRING,
    description="Load example dataset (1 or 2)",
    choices=["1", "2", None],
))

DATABASE_PARAM = register_param(QueryParam(
    key="database",
    default="metabolomicspanrepo_index_latest",
    param_type=ParamType.STRING,
    description="MASST database to search",
))

ANALOG_PARAM = register_param(QueryParam(
    key="analog",
    default=False,
    param_type=ParamType.BOOLEAN,
    description="Enable analog search",
))

PRECURSOR_TOL_PARAM = register_param(QueryParam(
    key="precursor_tol",
    default=0.02,
    param_type=ParamType.FLOAT,
    description="Precursor m/z tolerance (Da)",
    validator=lambda v: 0.001 <= v <= 1.0,
))

FRAGMENT_TOL_PARAM = register_param(QueryParam(
    key="fragment_tol",
    default=0.02,
    param_type=ParamType.FLOAT,
    description="Fragment m/z tolerance (Da)",
    validator=lambda v: 0.001 <= v <= 1.0,
))

MIN_COS_PARAM = register_param(QueryParam(
    key="min_cos",
    default=0.7,
    param_type=ParamType.FLOAT,
    description="Minimum cosine similarity threshold",
    validator=lambda v: 0.0 <= v <= 1.0,
))

ORGANISM_PARAM = register_param(QueryParam(
    key="organism",
    default="Humans",
    param_type=ParamType.STRING,
    description="Organism filter for results",
    choices=["Humans", "Rodents", "All organisms", "Manual Entry"],
))


# =============================================================================
# API FUNCTIONS
# =============================================================================

def load_all_params() -> dict[str, Any]:
    """
    Load all registered parameters from URL query params.
    
    Returns:
        Dict mapping parameter keys to decoded values
    """
    url_params = st.query_params.to_dict()
    result = {}
    
    for key, param in PARAM_REGISTRY.items():
        if key in url_params:
            result[key] = param.decode(url_params[key])
        else:
            result[key] = param.default
    
    return result


def sync_params(values: dict[str, Any], clear_first: bool = False):
    """
    Sync values to URL query parameters.
    
    Args:
        values: Dict of {param_key: value} to sync
        clear_first: Clear all params before syncing
    """
    if clear_first:
        st.query_params.clear()
    
    for key, value in values.items():
        if key not in PARAM_REGISTRY:
            continue
        
        param = PARAM_REGISTRY[key]
        encoded = param.encode(value)
        
        if encoded is not None and encoded != "":
            st.query_params[key] = encoded
        elif key in st.query_params:
            del st.query_params[key]
    
    # Special handling: if example is set, remove usi from URL
    if 'example' in values and values['example'] in ('1', '2'):
        if 'usi' in st.query_params:
            del st.query_params['usi']


def sync_param(key: str, value: Any):
    """Sync a single parameter to URL."""
    sync_params({key: value})


def get_param_schema() -> list[dict]:
    """
    Export schema for all registered parameters.
    
    Useful for external tools to discover supported parameters.
    
    Returns:
        List of parameter schemas (OpenAPI-style)
    """
    return [param.to_schema() for param in PARAM_REGISTRY.values()]


def get_defaults() -> dict[str, Any]:
    """Get default values for all parameters."""
    return {key: param.default for key, param in PARAM_REGISTRY.items()}


def get_default(key: str) -> Any:
    """Get default value for a specific parameter."""
    if key in PARAM_REGISTRY:
        return PARAM_REGISTRY[key].default
    return None


def clear_all_params():
    """Clear all query parameters from URL."""
    st.query_params.clear()
