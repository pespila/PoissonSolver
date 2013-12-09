#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_PAARE 255

void print_location(char *);
char *getdata();
char *Strdup(const char *);
void hex2ascii(char *);
char convert(char *);
struct CGI_DATEN *erstellen(char *);
void printf_error(char *);

struct CGI_DATEN {
   char *variable;
   char *wert;
   struct CGI_DATEN *next;
};

struct CGI_DATEN *ende = NULL;

/* Weiterleitung zu einer URL;
 * url ist die URL, wohin Sie den User weiterleiten.
 */
void print_location(char *url) {
   printf("Location: %s\n", url);
   /* für den Fall, dass ein alter Browser keine
      automatische Weiterleitung unterstützt */
   printf("Content-Type: text/html\n\n");
   printf("<html><head>\n");
   printf("<title>Weiterleitung zu %s</title>\n",url);
   printf("</head><body>\n");
   printf("Weiter gehts <a href=\"%s\">hier</a>",url);
   printf("</body></html>\n");
}

/*
 *  Funktion liest Daten in der POST- oder GET-Methode ein.
 *  Rückgabewert: String puffer mit den Daten
 *  bei Fehler  : NULL
 */
char *getdata(void) {
   unsigned long size;
   char *puffer = NULL;
   char *request = getenv("REQUEST_METHOD");
   char *cont_len;
   char *cgi_string;

   /* zuerst die Request-Methode überprüfen */
   if(  NULL == request )
      return NULL;
   else if( strcmp(request, "GET") == 0 ) {
      /* Die Methode GET -> Query-String abholen */
      cgi_string = getenv("QUERY_STRING");
      if( NULL == cgi_string )
         return NULL;
      else {
         puffer = (char *) Strdup(cgi_string);
         return puffer; /* Rückgabewert an den Aufrufer */
      }
   }
   else if( strcmp(request, "POST") == 0 ) {
      /* die Methode POST -> Länge des Strings
       * ermitteln (CONTENT_LENGTH) */
      cont_len = getenv("CONTENT_LENGTH");
      if( NULL == cont_len)
         return NULL;
      else {
         /* String CONTENT_LENGTH in unsigned long umwandeln */
         size = (unsigned long) atoi(cont_len);
         if(size <= 0)
            return NULL; /* Keine Eingabe!?!? */
      }
      /* jetzt lesen wir die Daten von stdin ein */
      puffer =(char *) malloc(size+1);
      if( NULL == puffer )
         return NULL;
      else {
         if( NULL == fgets(puffer, size+1, stdin) ) {
            free(puffer);
            return NULL;
         }
         else  /* Rückgabewerte an den Aufrufer */
            return puffer;
      }
   }

   /* Weder die GET-Methode noch die POST-Methode wurden verwendet. */
   else
      return NULL;
}

/*  Da die Funktion strdup() in der Headerdatei <string.h> keine
 *  ANSI-C-Funktion ist, schreiben wir eine eigene.
 */
char *Strdup(const char *str) {
   char *p;
   if(NULL == str)
      return NULL;
   else {
      p = (char *)malloc(strlen(str)+1);
      if(NULL == p)
         return NULL;
      else
         strcpy(p, str);
   }
   return p;
}

/* Wandelt einzelne Hexzeichen (%xx) in ASCII-Zeichen
 * und kodierte Leerzeichen (+) in echte Leerzeichen um. */
void hex2ascii(char *str) {
   int x,y;

   for(x=0,y=0; str[y] != '\0'; ++x,++y) {
      str[x] = str[y];
      /* Ein hexadezimales Zeichen? */
      if(str[x] == '%')  {
         str[x] = convert(&str[y+1]);
         y += 2;
      }
      /* Ein Leerzeichen? */
      else if( str[x] == '+')
         str[x]=' ';
   }
   /* geparsten String sauber terminieren */
   str[x] = '\0';
}




/* Funktion konvertiert einen String von zwei hexadezimalen
 * Zeichen und gibt das einzelne dafür stehende Zeichen zurück.
 */
char convert(char *hex) {
   char ascii;

   /* erster Hexawert */
   ascii =
   (hex[0] >= 'A' ? ((hex[0] & 0xdf) - 'A')+10 : (hex[0] - '0'));
   ascii <<= 4; /* Bitverschiebung schneller als ascii*=16 */
   /* zweiter Hexawert */
   ascii +=
   (hex[1] >= 'A' ? ((hex[1] & 0xdf) - 'A')+10 : (hex[1] - '0'));
   return ascii;
}

/* Liste aus Variable/Wert-Paaren erstellen
 * Rückgabewert: Anfangsadresse der Liste
 * Bei Fehler: NULL
 */
struct CGI_DATEN *erstellen(char *str) {
   char* s;
   char* res;
   /* Irgendwo gibt es auch eine Grenze, hier sind
      MAX_PAARE erlaubt. */
   char *paare[MAX_PAARE];
   struct CGI_DATEN *ptr_daten = NULL;
   struct CGI_DATEN *ptr_anfang = NULL;
   int i=0, j=0;

   /* Zuerst werden die Variablen/Werte-Paare anhand des Zeichens
    * '&' getrennt, sofern es mehrere sind. */
    s=str;
    res=strtok(s,"&");
    while( res != NULL && i < MAX_PAARE) {
       /* Wert von res dynamisch in char **pair speichern */
       paare[i] = (char *)malloc(strlen(res)+1);
       if(paare[i] == NULL)
          return NULL;
       paare[i] = res;
       res=strtok(NULL,"&");
       i++;
    }


   /* Jetzt werden die Variablen von den Werten getrennt und
    * an die Struktur CGI_DATEN übergeben. */
   while ( i > j )  { /* Das erste Element? */
      if(ptr_anfang == NULL) {
         ptr_anfang =(struct CGI_DATEN *)
           malloc(sizeof (struct CGI_DATEN *));
         if( ptr_anfang == NULL )
            return NULL;
         res = strtok( paare[j], "=");
         if(res == NULL)
            return NULL;
         ptr_anfang->variable = (char *)
           malloc(strlen(res)+1);
         if( ptr_anfang->variable == NULL )
            return NULL;
         ptr_anfang->variable = res;
         res = strtok(NULL, "\0");
         if(res == NULL)
            return NULL;
         ptr_anfang->wert = (char *) malloc(strlen(res)+1);
         if( ptr_anfang->wert == NULL )
            return NULL;
         ptr_anfang->wert = res;
         /* printf("%s %s<br>",
          * ptr_anfang->variable, ptr_anfang->wert); */
         ptr_anfang->next =(struct CGI_DATEN *)
           malloc(sizeof (struct CGI_DATEN *));
         if(ptr_anfang->next == NULL)
            return NULL;
         ptr_daten = ptr_anfang->next;
         j++;
      }
      else { /* die restlichen Elemente */
         res = strtok( paare[j], "=");
         if(res == NULL)
            return NULL;
         ptr_daten->variable =(char *)
           malloc(strlen(res)+1);
         if(ptr_daten->variable == NULL)
            return NULL;
         ptr_daten->variable = res;
         res = strtok(NULL, "\0");
         if(res == NULL)
            return NULL;
         ptr_daten->wert =(char *) malloc(strlen(res)+1);
         if(ptr_daten->wert == NULL)
            return NULL;
         ptr_daten->wert = res;
         /* printf("%s %s<br>",
          * ptr_daten->variable,  ptr_daten->wert); */
         ptr_daten->next = (struct CGI_DATEN *)
           malloc(sizeof (struct CGI_DATEN *));
         if( ptr_daten->next == NULL )
            return NULL;
         ptr_daten = ptr_daten->next;
         j++;
      }
   }
   ende = ptr_daten;
   /* Anfangsadresse der Liste struct CGI_DATEN zurückgeben */
   return ptr_anfang;
}

void loeschen(struct CGI_DATEN *daten) {
   struct CGI_DATEN *next = NULL;

   while(daten != ende) {
      next = daten->next;
      if(daten->variable != NULL)
         free(daten);
      daten=next;
   }
}

void printf_error(char *str) {
   printf("Content-Type: text/html\n\n");
   printf("<html><head>\n");
   printf("<title>CGI-Fehlermeldung</title>\n");
   printf("</head><body>\n");
   printf("%s",str);
   printf("</body></html>\n");
}

int main(void) {
   char *str;
   struct CGI_DATEN *cgi;
   struct CGI_DATEN *free_cgi;
   FILE *f;


   /* Eingabe einlesen */
   str = getdata();
   if(str == NULL) {
      printf_error("Fehler beim Einlesen von der "
                   "Formulareingabe");
      return EXIT_FAILURE;
   }
   /* Hexzeichen in ASCII-Zeichen konvertieren und aus '+'
    * Leerzeichen machen */
   hex2ascii(str);
   /* Liste der Formualardaten erstellen */
   cgi = erstellen(str);
   free_cgi = cgi;
   if (cgi == NULL) {
      printf_error("Fehler beim Erstellen der "
                   "Variablen/Werte-Liste\n");
      return EXIT_FAILURE;
   }

   /* Datei zum Schreiben öffnen */
   /* Bitte den Pfad anpassen: beispielsweise unter SUSE Linux:
    * f = fopen("/srv/www/htdocs/gaeste.html", "r+");
    * und WICHTIG: Schreibrechte auf diese Datei vergeben
    */
   f = fopen("gaeste.html", "r+");
   if(f == NULL) {
      printf_error("Konnte Datei gaeste.html nicht zum "
                   "Schreiben oeffnen\n");
      return EXIT_FAILURE;
   }
   else {
      /* Stream vor </body></html> */
      fseek(f, -14, SEEK_END);
      fprintf(f, "<hr><br>"); /* Eine horizontale Linie */
      /* Name */
      if(cgi->wert != NULL)
         fprintf(f, "Name: %s E-Mail: ",cgi->wert);
      cgi = cgi->next;
      /* Mailadresse */
      if(cgi->wert != NULL)
         fprintf(f, "<a href=\"mailto:%s\">%s</a> ",
           cgi->wert,cgi->wert);
      cgi = cgi->next;
      /* Bewertung */
      if(cgi->wert != NULL)
         fprintf(f, "Bewertung : %s",cgi->wert);
      cgi = cgi->next;
      /* Eintrag */
      if( cgi->wert != NULL) {
         fprintf(f, "<p><b>Der Eintrag : </b>");
         fprintf(f, "%s",cgi->wert);
      }
      cgi = cgi->next;
      /* Programmierkenntnis(se) */
      if(cgi->wert != NULL) {
         fprintf(f, "<br><br>Programmierkenntnisse : ");
         while(cgi->wert != NULL  &&
           strcmp(cgi->variable,"programmieren") == 0 ) {
            fprintf(f, "%s ",cgi->wert);
            cgi = cgi->next;
         }
      }
      fprintf(f, "</p></body></html>");
      fclose(f);
   }
   /* Speicher wieder freigeben */
   loeschen(free_cgi);
   /* Auch hier müssen Sie die Pfadangabe ggf. anpassen. */
   print_location("http://localhost/gaeste.html");
   return EXIT_SUCCESS;
}