import os
from sqlalchemy import create_engine

user = os.getenv('MYSQL_USER')
password = os.getenv('MYSQL_PASSWORD')
host = os.getenv('MYSQL_HOST')
port = os.getenv('MYSQL_PORT')
database = os.getenv('MYSQL_DATABASE')

def create_sql_engine():
  """Create a SQLAlchemy engine for the IEDB MySQL backend."""
  return create_engine(f"mysql+mysqlconnector://{user}:{password}@{host}:{port}/{database}")