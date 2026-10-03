# MCL Explorer frontend

Next.js research interface for Master Cancer Landscape.

## Local run

```powershell
cd C:\Users\Bogdan\Desktop\master-cancer-landscape\apps\explorer
npm install
$env:NEXT_PUBLIC_MCL_API_URL="http://127.0.0.1:8000"
npm run dev
```

Open `http://localhost:3000`.

The FastAPI service must be running on port 8000.
