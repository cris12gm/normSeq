from django.contrib import admin
from django.urls import include, re_path
from django.conf import settings
from django.conf.urls.static import static
from .views import Query

urlpatterns = [
    re_path('', Query.as_view(), name='query'),
    ]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)