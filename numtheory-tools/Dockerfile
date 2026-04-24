# syntax=docker/dockerfile:1

# Маленький и быстрый базовый образ
FROM python:3.12-slim

# Рабочая директория внутри контейнера
WORKDIR /app

# Копируем сначала requirements (кешируется, если зависимости не меняются)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Копируем сам код
COPY . .

# Фронтовой порт
EXPOSE 8080

# Запуск
# host=0.0.0.0 — чтобы сервер был доступен извне контейнера
CMD ["python", "app.py"]
